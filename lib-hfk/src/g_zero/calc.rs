use std::str::FromStr;

use log::*;
use num_traits::{One, Zero};
use yui_core::lc::{LcKey, Lc};
use yui_core::num::FF2;
use yui_core::poly::MultiVar;
use yui_homology::{ToTableString, Summand};
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};
use yui_matrix::sparse::pluq::{solve_pluq, solve_pluq_incr};

use crate::{GCGen, GridComplex, GridDiagram, HFK};
use crate::misc::lanczos::{find_ng_certificate, find_solution};
use crate::grid_complex::{GCFiltered, GCInfGen};
use crate::grid_complex::chain::GCInfChainExt;
use super::{Candidates, Region};

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct G0CalcConfig {
    pub filter_regions: bool,
    pub filter_gens: bool,
    pub accept_if_unique: bool,
    pub try_fast_ng: bool,
    pub mode: Mode,
    pub max_size: usize,
    pub max_piv: usize,
    pub chunk: usize,
}

impl Default for G0CalcConfig {
    fn default() -> Self {
        Self {
            filter_regions: true,
            filter_gens: true,
            accept_if_unique: false,
            try_fast_ng: false,
            mode: Mode::default(),
            max_size: 1_000_000,
            max_piv: 500_000,
            chunk: 100_000,
        }
    }
}

/// The solver for the linear systems over `𝔽₂`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Mode {
    /// Block Lanczos: small memory, randomized but with checked answers.
    #[default]
    Lanczos,
    /// (incremental) PLUQ: deterministic, but may fill in on large systems.
    Pluq,
}

impl FromStr for Mode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "lanczos" => Ok(Mode::Lanczos),
            "pluq" => Ok(Mode::Pluq),
            _ => Err(format!("unknown mode `{s}` (expected `lanczos` or `pluq`)")),
        }
    }
}

// note: one Lanczos run almost surely finds a certificate or a solution when one exists.
const LANCZOS_ROUNDS: usize = 1;

/// Returned by `find_homologous` when the system exceeds `max_size` and
/// neither the trivial cycle-membership check nor the U-specialization
/// fast-path concluded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum G0CalcError { Resigned }

pub struct G0Calc {
    gc: GridComplex<GCFiltered>,
    z: Lc<GCInfGen, FF2>,
    candidates: Candidates,
    config: G0CalcConfig
}

impl G0Calc {
    /// `g`, `g4` and `tau` are known values, if any: `g` is checked against HFK,
    /// and `g4`, `tau` are used to filter the candidates.
    pub fn init(gd: &GridDiagram, g: Option<usize>, g4: Option<usize>, tau: Option<isize>, config: G0CalcConfig) -> Self {
        info!("compute HFK ...");
        let hfk = HFK::from(gd);
        info!("HFK:\n{}", hfk.to_table_string());
        info!("g = {}", hfk.genus());

        Self::validate(&hfk, g, g4, tau);

        let candidates = Self::make_candidates(&hfk, g4, tau, config);
        Self::init_with_candidates(gd, candidates, config)
    }

    /// Tests only the given candidates, in the given order.
    pub fn init_with_candidates(gd: &GridDiagram, candidates: Candidates, config: G0CalcConfig) -> Self {
        info!("prepare generators ...");
        let gc = Self::build_complex(gd, &candidates, config.filter_gens);
        info!("  #generators: {}", gc.gens.len());

        info!("compute homological generator..");
        let z = gc.hom_gen();
        
        if z.is_zero() { 
            info!("  z = 0.");
        } else {
            info!("  support of z: \n{}", z.support_grid().to_table_string());
        }

        G0Calc { gc, z, candidates, config }
    }

    fn validate(hfk: &HFK, g: Option<usize>, g4: Option<usize>, tau: Option<isize>) {
        let genus = hfk.genus() as usize;

        if let Some(g) = g {
            assert_eq!(g, genus, "expected genus g = {g}, but computed g = {genus}.");
        }

        // MEMO: replace with if-let chain
        if let Some(g4) = g4 && let Some(tau) = tau { 
            let t_abs = tau.unsigned_abs();
            assert!(t_abs <= g4, "Invalid inequality: |tau| = {t_abs} > g4 = {g4}");
        }
        
        if let Some(g4) = g4 { 
            assert!(g4 <= genus, "Invalid inequality: g4 = {g4} > g = {genus}");
        }

        if let Some(tau) = tau { 
            assert!(hfk[(0, tau)] > 0, "tau = {tau} is incompatible with HFK.");
        }
    }

    fn make_candidates(hfk: &HFK, g4: Option<usize>, tau: Option<isize>, config: G0CalcConfig) -> Candidates { 
        let g = hfk.genus();
        let mut candidates = Candidates::generate_from(g as usize);

        debug!("initial candidate regions: {}", candidates.regions.len());

        if config.filter_regions { 
            candidates.desymmetrize();
            debug!("  desymmetrize -> {}", candidates.regions.len());

            candidates.filter_by_hfk(hfk);
            debug!("  filter by hfk -> {}", candidates.regions.len());

            if let Some(g4) = g4 { 
                candidates.filter_by_g4(g4);
                debug!("  filter by g4={g4} -> {}", candidates.regions.len());
            }

            if let Some(tau) = tau { 
                candidates.filter_by_tau(tau);
                debug!("  filter by tau={tau} -> {}", candidates.regions.len());
            }
        }

        candidates.sort();

        info!("candidate regions: {}", candidates.regions.len());
        for r in &candidates.regions {
            info!("  {r}");
        }

        if config.accept_if_unique && candidates.regions.len() == 1 { 
            let r = candidates.regions.pop().unwrap();
            info!("accept unique region:");
            info!("  {r}");

            candidates.accept(&r);
        }

        candidates
    }

    fn build_complex(gd: &GridDiagram, candidates: &Candidates, filter_gens: bool) -> GridComplex<GCFiltered> {
        let gd = gd.clone();
        if !filter_gens { return GridComplex::new(gd); }

        GridComplex::new_filtered(gd, |x| 
            candidates.regions.iter().any(|r| { 
                let s = r.right_end(); // r is to be shifted by (-s, -s).
                let m = x.maslov_deg;
                if m < -2 * s { return false; }

                // Consider U^a x of total degree e such that
                //
                //    M(U^a x) = -2s or -2s + 1.
                //
                // Then A(U^a x) = -e + A(x), while alg(U^a x) = -a₀ ranges
                // over -e ..= 0 as a₀ varies. Since r is down-closed, the
                // inflation most likely to lie outside is the one with a₀ = 0,
                // so it suffices to test alg = 0.

                let e = (m + 2 * s) / 2; // >= 0
                !r.contains((s, s - e + x.alex_deg)) // <=> !r.shift(-s, -s).contains(0, -e + x.alex_deg)
            })
        )
    }

    pub fn complex(&self) -> &GridComplex<GCFiltered> {
        &self.gc
    }

    /// Processes all remaining candidates and populates `candidates.accepted`.
    pub fn run(&mut self) {
        while !self.candidates.regions.is_empty() {
            let r = self.candidates.regions.remove(0);
            self.process(&r);
        }
    }

    pub fn accepted(&self) -> &[Region] {
        &self.candidates.accepted
    }

    pub fn resigned(&self) -> &[Region] {
        &self.candidates.resigned
    }

    // Tests whether `r` is realizable: tries to find a cycle homologous to z
    // inside r. Shifts both r and z by s = right_end(r) so that the shifted
    // region has right_end = 0, which is required by find_homologous.
    fn process(&mut self, r: &Region) {
        info!("try: {r} (remain: {})", self.candidates.regions.len());

        if r.contains_chain(&self.z) {
            debug!("  z contained in {r}.");
            info!("OK: {r}");
            self.candidates.accept(r);
            return
        }

        let s = r.right_end().max(0);
        let r_shifted = r.shift(-s, -s);
        let z_shifted = self.shift_z(s as usize);

        match self.find_homologous(&z_shifted, &r_shifted) {
            Ok(Some(_)) => {
                info!("OK: {r}");
                self.candidates.accept(r);
            }
            Ok(None) => {
                info!("NG: {r}");
                self.candidates.reject(r);
            }
            Err(G0CalcError::Resigned) => {
                info!("RESIGN: {r}");
                self.candidates.resign(r);
            }
        }
    }

    // Returns z multiplied by U_0^s, shifting all generator coordinates by
    // (-s, -s) to match r.shift(-s, -s).
    fn shift_z(&self, s: usize) -> Lc<GCInfGen, FF2> {
        let u = MultiVar::from((0, s));
        self.z.mul(&u)
    }

    /// Tests whether `z` has a homologous cycle contained in the region `r`.
    ///
    /// Returns:
    /// - `Ok(Some(z'))` if a homologous cycle was found.
    /// - `Ok(None)` if the system is unsolvable (definitive NG).
    /// - `Err(Resigned)` if the system exceeds `max_size` and no fast-path
    ///   could decide it.
    fn find_homologous(
        &self,
        z: &Lc<GCInfGen, FF2>,
        r: &Region,
    ) -> Result<Option<Lc<GCInfGen, FF2>>, G0CalcError> {
        assert!(z.is_maslov_homogeneous());

        let k = z.maslov_degree();
        let z0 = r.cutoff(z);

        if z0.is_zero() { 
            return Ok(Some(Lc::zero()))
        }

        let pts_w = self.gc.complementary_points(k + 1, r);
        let pts_z = self.gc.complementary_points(k, r);
        let n_w = self.gc.count_inf_gens_in(k + 1, &pts_w);
        let n_z = self.gc.count_inf_gens_in(k, &pts_z);

        debug!("linear system size: {:?}", (n_z, n_w));

        // try the NG fast-path
        if self.config.try_fast_ng {
            let n_raw_w = self.gc.count_raw_gens_in(k + 1, &pts_w);
            let n_raw_z = self.gc.count_raw_gens_in(k, &pts_z);
            let smaller = n_raw_w < n_w || n_raw_z < n_z;
            let tractable = n_raw_w < self.config.max_size && n_raw_z < self.config.max_size;
            if smaller && tractable && !self.pass_specialized(k, &pts_w, &pts_z, &z0, r) {
                return Ok(None);
            }
        }

        // Resign if the inflated system is too large for the full solve.
        if n_w > self.config.max_size || n_z > self.config.max_size {
            debug!("resigning: size {} > max_size {}", n_w.max(n_z), self.config.max_size);
            return Err(G0CalcError::Resigned);
        }

        debug!("create inf gens ...");
        let gen_w = self.gc.inf_gens_in(k + 1, &pts_w); // space where w lives
        let gen_z = self.gc.inf_gens_in(k, &pts_z);     // space where z lives mod r

        debug!("create matrix ...");
        let (a, y) = create_ay(
            &gen_w, &gen_z, &z0,
            |x| r.cutoff(&self.gc.differentiate(x)),
        );
        debug!("  a: {:?} (nnz: {})", a.shape(), a.nnz());

        let Some(x_sol) = self.find_sol(&a, &y)? else {
            return Ok(None);
        };
        let w = gen_w.devectorize(&x_sol);

        // z - d(w) is homologous to z and lies inside r.
        let mut z1 = z.clone();
        z1 -= &self.gc.differentiate(&w);
        Ok(Some(z1))
    }

    /// `Ok(Some(x))`: a solution of `a·x = y`, `Ok(None)`: proven unsolvable,
    /// `Err(Resigned)`: undecided (Lanczos found neither a certificate nor a solution).
    fn find_sol(&self, a: &SpMat<FF2>, y: &SpVec<FF2>) -> Result<Option<SpVec<FF2>>, G0CalcError> {
        debug!("solve linear system ...");

        let x = match self.config.mode {
            Mode::Lanczos => {
                if find_ng_certificate(a, y, LANCZOS_ROUNDS).is_some() {
                    return Ok(None);
                }
                find_solution(a, y, LANCZOS_ROUNDS).ok_or(G0CalcError::Resigned)?
            }
            Mode::Pluq => {
                let x = if a.n_rows() < self.config.chunk {
                    solve_pluq(a, y)
                } else {
                    solve_pluq_incr(a, y, self.config.max_piv, self.config.chunk)
                };
                let Some(x) = x else {
                    return Ok(None);
                };
                x
            }
        };

        debug!("validate solution..");
        assert_eq!(&(a * &x), y);
        debug!("  ok.");

        Ok(Some(x))
    }

    fn pass_specialized(
        &self,
        k: isize,
        pts_w: &[(isize, isize)],
        pts_z: &[(isize, isize)],
        z0: &Lc<GCInfGen, FF2>,
        r: &Region,
    ) -> bool {
        debug!("try NG fast-path ...");

        let n = self.gc.grid_number() as usize;
        let raw_gen_w = self.gc.raw_gens_in(k + 1, pts_w);
        let raw_gen_z = self.gc.raw_gens_in(k, pts_z);

        if self.no_sol_specialized(&raw_gen_w, &raw_gen_z, z0, r, &vec![FF2::one(); n]) {
            debug!("  NG via U = 1");
            return false;
        }

        debug!("passed NG fast-path.");
        true
    }

    fn no_sol_specialized(
        &self,
        raw_gen_w: &Summand<GCGen, FF2>,
        raw_gen_z: &Summand<GCGen, FF2>,
        z0: &Lc<GCInfGen, FF2>,
        r: &Region,
        c: &[FF2],
    ) -> bool {
        let lift = |x: &GCGen| Lc::from(GCInfGen::from(x.clone()));
        let proj = move |g: &GCInfGen| if g.u.eval::<FF2>(c).is_zero() {
            Lc::zero()
        } else {
            Lc::from(g.x.clone())
        };

        debug!("create specialized matrix ...");
        let (a_c, y_c) = create_ay(
            raw_gen_w, raw_gen_z, &z0.apply(proj),
            |x| r.cutoff(&self.gc.differentiate(&x.apply(lift))).apply(proj),
        );

        debug!("  specialized: {:?} (nnz = {})", a_c.shape(), a_c.nnz());

        // note: only "proven unsolvable" matters here, so Lanczos skips the search for a solution.
        match self.config.mode {
            Mode::Lanczos => find_ng_certificate(&a_c, &y_c, LANCZOS_ROUNDS).is_some(),
            Mode::Pluq => matches!(self.find_sol(&a_c, &y_c), Ok(None)),
        }
    }
}

fn create_ay<X>(
    gen_w: &Summand<X, FF2>,
    gen_z: &Summand<X, FF2>,
    y_target: &Lc<X, FF2>,
    f: impl Fn(&Lc<X, FF2>) -> Lc<X, FF2>,
) -> (SpMat<FF2>, SpVec<FF2>)
where
    X: LcKey + Into<GCInfGen>,
{
    let a = gen_w.make_matrix(gen_z, |x| f(x));
    let y = gen_z.vectorize(y_target);
    (a, y)
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;

    use super::*;
    use crate::GridDiagram;
    use crate::g_zero::Corner;

    fn config() -> G0CalcConfig {
        G0CalcConfig { 
            filter_regions: false, 
            filter_gens: false, 
             ..Default::default() 
        }
    }

    #[test]
    fn test_find_homologous_unknot_already_inside() {
        // For the unknot the hom_gen lies at (alg=0, alex=0) ∈ R(0,0),
        // so z0 = cutoff(z, R(0,0)) = 0.  The system is trivially satisfied
        // with w = 0, and the function returns z unchanged.
        let gd = GridDiagram::unknot();
        let cfg = config();
        let calc = G0Calc::init(&gd, None, None, None, cfg);

        let z = &calc.z;
        assert!(!z.is_zero());

        let r = Region::simple(0, 0);
        assert!(r.contains_chain(z));
        
        let z1 = calc.find_homologous(z, &r).unwrap().expect("expected Some");
        assert!(z1.is_zero());
    }

    #[test]
    fn test_find_homologous_unknot_none() {
        // R[(-1,0), (0,-1)] lies strictly below the hom_gen at (0,0).
        let gd = GridDiagram::unknot();
        let cfg = config();
        let calc = G0Calc::init(&gd, None, None, None, cfg);
        let z = &calc.z;

        assert!(!z.is_zero());

        let r = Region::new([Corner::new(-1, 0), Corner::new(0, -1)]);
        assert_eq!(calc.find_homologous(z, &r), Ok(None));
    }

    #[test]
    fn test_find_homologous_trefoil() {
        // For the trefoil (genus 1), R(0,1) is the expected realizable region.
        // find_homologous should return Some(z1) where z1 is a cycle inside R(0,1).
        let gd = GridDiagram::trefoil();
        let cfg = config();
        let calc = G0Calc::init(&gd, None, None, None, cfg);
        let z = &calc.z;
        assert!(!z.is_zero());

        let r = Region::simple(0, 1);
        assert!(r.contains_chain(z));

        let z1 = calc.find_homologous(z, &r).unwrap().expect("expected Some for trefoil R(0,1)");
        assert!(z1.is_zero());
    }

    fn accepted(gd: GridDiagram) -> Vec<Region> {
        let cfg = config();
        let mut calc = G0Calc::init(&gd, None, None, None, cfg);
        calc.run();

        let mut rs = calc.accepted().to_vec();
        rs.sort_by_key(|r| r.right_end());
        rs
    }

    #[test]
    fn test_unknot() {
        let rs = accepted(GridDiagram::test_data("0_1").unwrap());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_trefoil() {
        let rs = accepted(GridDiagram::test_data("3_1").unwrap());
        assert_eq!(rs.len(), 2);
        assert!(rs.contains(&Region::simple(0, 1)));
        assert!(rs.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_trefoil_mirror() {
        let rs = accepted(GridDiagram::test_data("3_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([Corner::new(-1, 0), Corner::new(0, -1)])]);
    }

    #[test]
    fn test_figure8() {
        let rs = accepted(GridDiagram::test_data("4_1").unwrap());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_figure8_mirror() {
        let rs = accepted(GridDiagram::test_data("4_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_5_1() {
        let rs = accepted(GridDiagram::test_data("5_1").unwrap());
        assert_eq!(rs.len(), 3);
        assert!(rs.contains(&Region::simple(0, 2)));
        assert!(rs.contains(&Region::simple(1, 1)));
        assert!(rs.contains(&Region::simple(2, 0)));
    }

    #[test]
    fn test_5_1_mirror() {
        let rs = accepted(GridDiagram::test_data("5_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([
            Corner::new(-2, 0),
            Corner::new(-1, -1),
            Corner::new(0, -2),
        ])]);
    }

    #[test]
    fn test_5_2() {
        let rs = accepted(GridDiagram::test_data("5_2").unwrap());
        assert_eq!(rs.len(), 2);
        assert!(rs.contains(&Region::simple(0, 1)));
        assert!(rs.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_5_2_mirror() {
        let rs = accepted(GridDiagram::test_data("5_2").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([Corner::new(-1, 0), Corner::new(0, -1)])]);
    }

    #[test]
    fn test_6_1() {
        let rs = accepted(GridDiagram::test_data("6_1").unwrap());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_6_2() {
        let rs = accepted(GridDiagram::test_data("6_2").unwrap());
        assert_eq!(rs.len(), 2);
        assert!(rs.contains(&Region::simple(0, 1)));
        assert!(rs.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_6_3() {
        let rs = accepted(GridDiagram::test_data("6_3").unwrap());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    // G0(K) is closed under the reflection (i, j) ↦ (j, i), which follows from
    // the symmetry CFK^∞(K) ≃ (CFK^∞(K))^r. See [SS2019] §5.1.
    #[test]
    fn test_reflection_symmetry() {
        for name in ["3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3"] {
            let gd = GridDiagram::test_data(name).unwrap();
            let rs = accepted(gd.clone());

            // every accepted region must be non-empty and reflection-closed as a set
            assert!(!rs.is_empty(), "{name}: no accepted region");
            for r in &rs {
                assert!(rs.contains(&r.transposed()), "{name}: {r} has no reflection");
            }
        }
    }
}
