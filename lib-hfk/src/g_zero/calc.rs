use log::*;
use num_traits::Zero;
use yui_core::lc::{LcKey, Lc};
use yui_core::num::FF2;
use yui_core::poly::MultiVar;
use yui_homology::{ToTableString, Summand};
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};
use yui_matrix::sparse::pluq::{solve_pluq, solve_pluq_incr};

use crate::{GridComplex, GridDiagram, HFK};
use crate::misc::lanczos::{find_ng_certificate, find_solution};
use crate::grid_complex::{GCFiltered, GCInfGen};
use crate::grid_complex::chain::GCInfChainExt;
use super::{Candidates, Region};
use super::complex_ext::project;

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct G0CalcConfig {
    pub mode: Mode,
    pub filter_regions: bool,
    pub filter_gens: bool,
    /// Accept the last candidate without a test when nothing is accepted or resigned so far.
    /// Valid only when the candidates contain `G0(K)`, as those of `G0Calc::init` do.
    pub accept_if_unique: bool,
    pub try_fast_ng: bool,
    pub max_size: usize,
}

impl Default for G0CalcConfig {
    fn default() -> Self {
        Self {
            mode: Mode::default(),
            filter_regions: true,
            filter_gens: true,
            accept_if_unique: false,
            try_fast_ng: false,
            max_size: 1_000_000,
        }
    }
}

/// The solver for the linear systems over `𝔽₂`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// Block Lanczos: small memory, randomized but with checked answers.
    /// Tries up to `rounds` random starts.
    Lanczos { rounds: usize },
    /// PLUQ: deterministic, but may fill in on large systems.
    /// Large systems are solved incrementally, in steps of `chunk`.
    Pluq { chunk: usize },
}

impl Mode {
    pub fn lanczos(rounds: usize) -> Self {
        Mode::Lanczos { rounds }
    }

    pub fn pluq(chunk: usize) -> Self {
        Mode::Pluq { chunk }
    }
}

impl Default for Mode {
    // note: one Lanczos run almost surely finds a certificate or a solution when one exists.
    fn default() -> Self {
        Mode::lanczos(1)
    }
}

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

    /// Regions expected to be realizable, as a hint for the order of the tests: the candidates
    /// among them (up to transposition) are tested first, searching for a solution first and
    /// without the NG fast-path. The results do not depend on the hint.
    pub fn set_expected(&mut self, regions: Vec<Region>) {
        self.candidates.set_expected(regions);
    }

    /// Processes all remaining candidates and populates `candidates.accepted`: the expected
    /// ones first, searching for a solution first, and then the rest as usual.
    pub fn run(&mut self) {
        while let Some((r, expected)) = self.candidates.take_next() {
            self.process(&r, expected);
        }
    }

    pub fn accepted(&self) -> &[Region] {
        &self.candidates.accepted
    }

    pub fn rejected(&self) -> &[Region] {
        &self.candidates.rejected
    }

    pub fn resigned(&self) -> &[Region] {
        &self.candidates.resigned
    }

    /// `G0(K)`, the accepted regions, unless a resigned region may still be minimal: one that
    /// contains no accepted region and lies in no rejected one.
    pub fn g0(&self) -> Result<Vec<Region>, G0CalcError> {
        let undecided = self.resigned().iter().any(|r| {
            !self.accepted().iter().any(|a| r.includes(a)) && !self.rejected().iter().any(|s| s.includes(r))
        });
        if undecided {
            Err(G0CalcError::Resigned)
        } else {
            Ok(self.accepted().to_vec())
        }
    }

    // Tests whether `r` is realizable: tries to find a cycle homologous to z
    // inside r. Shifts both r and z by s = right_end(r) so that the shifted
    // region has right_end = 0, which is required by find_homologous.
    fn process(&mut self, r: &Region, positive: bool) {
        info!("try: {r} (remain: {})", self.candidates.regions.len());

        // note: G0 is non-empty and among the candidates, and all the others are rejected.
        let c = &self.candidates;
        if self.config.accept_if_unique && c.regions.is_empty() && c.accepted.is_empty() && c.resigned.is_empty() {
            info!("OK: {r} (the last candidate)");
            self.candidates.accept(r);
            return
        }

        if r.contains_chain(&self.z) {
            debug!("  z contained in {r}.");
            info!("OK: {r}");
            self.candidates.accept(r);
            return
        }

        let s = r.right_end().max(0);
        let r_shifted = r.shift(-s, -s);
        let z_shifted = self.shift_z(s as usize);

        if positive {
            debug!("  {r} is expected: positive first.");
        }

        match self.find_homologous(&z_shifted, &r_shifted, positive) {
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

    /// Tests whether `z` has a homologous cycle contained in the region `r`. With `positive`,
    /// `r` is expected to be realizable: a solution is searched first, and the NG fast-path is skipped.
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
        positive: bool,
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

        // try the NG fast-path, unless r is expected to be realizable
        if self.config.try_fast_ng && !positive {
            let n_proj_w = self.gc.count_proj_gens_in(k + 1, &pts_w);
            let n_proj_z = self.gc.count_proj_gens_in(k, &pts_z);
            debug!("specialized system size: {:?}", (n_proj_z, n_proj_w));
            let smaller = n_proj_w < n_w || n_proj_z < n_z;
            let tractable = n_proj_w < self.config.max_size && n_proj_z < self.config.max_size;
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

        let Some(x_sol) = self.find_sol(&a, &y, positive)? else {
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
    /// With `positive`, Lanczos searches for a solution before a certificate.
    fn find_sol(&self, a: &SpMat<FF2>, y: &SpVec<FF2>, positive: bool) -> Result<Option<SpVec<FF2>>, G0CalcError> {
        debug!("solve linear system ...");

        let x = match self.config.mode {
            Mode::Lanczos { rounds } => {
                let certify = || find_ng_certificate(a, y, rounds).is_some();
                let solve = || find_solution(a, y, rounds);

                if positive {
                    match solve() {
                        Some(x) => x,
                        None if certify() => return Ok(None),
                        None => return Err(G0CalcError::Resigned),
                    }
                } else {
                    if certify() {
                        return Ok(None);
                    }
                    solve().ok_or(G0CalcError::Resigned)?
                }
            }
            Mode::Pluq { chunk } => {
                let x = if a.n_rows() < chunk {
                    solve_pluq(a, y)
                } else {
                    solve_pluq_incr(a, y, chunk, chunk)
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

        let gen_w = self.gc.proj_gens_in(k + 1, pts_w);
        let gen_z = self.gc.proj_gens_in(k, pts_z);

        if self.no_sol_specialized(&gen_w, &gen_z, z0, r) {
            debug!("  NG via Uᵢ = U₁ (i ≥ 1)");
            return false;
        }

        debug!("passed NG fast-path.");
        true
    }

    fn no_sol_specialized(
        &self,
        gen_w: &Summand<GCInfGen, FF2>,
        gen_z: &Summand<GCInfGen, FF2>,
        z0: &Lc<GCInfGen, FF2>,
        r: &Region,
    ) -> bool {
        let proj = |g: &GCInfGen| Lc::from(project(g));

        debug!("create specialized matrix ...");
        let (a_c, y_c) = create_ay(
            gen_w, gen_z, &z0.apply(proj),
            |x| r.cutoff(&self.gc.differentiate(x)).apply(proj),
        );

        debug!("  specialized: {:?} (nnz = {})", a_c.shape(), a_c.nnz());

        // note: only "proven unsolvable" matters here, so Lanczos skips the search for a solution.
        match self.config.mode {
            Mode::Lanczos { rounds } => find_ng_certificate(&a_c, &y_c, rounds).is_some(),
            Mode::Pluq { .. } => matches!(self.find_sol(&a_c, &y_c, false), Ok(None)),
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
    use std::collections::HashSet;

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

    // note: with accept_if_unique, the last candidate is accepted without a test. R(-1, 1) is
    // not realizable for the trefoil, but these candidates do not contain G0.
    #[test]
    fn test_accept_if_unique() {
        let gd = GridDiagram::trefoil();
        let accepted = |accept_if_unique| {
            let mut candidates = Candidates::new();
            candidates.add(Region::simple(0, 0));
            candidates.add(Region::simple(-1, 1));
            let mut calc = G0Calc::init_with_candidates(&gd, candidates, G0CalcConfig { accept_if_unique, ..config() });
            calc.run();
            calc.accepted().iter().cloned().collect::<HashSet<_>>()
        };
        assert!(accepted(false).is_empty());
        assert_eq!(accepted(true), HashSet::from([Region::simple(-1, 1), Region::simple(1, -1)]));
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
        
        let z1 = calc.find_homologous(z, &r, false).unwrap().expect("expected Some");
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
        assert_eq!(calc.find_homologous(z, &r, false), Ok(None));
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

        let z1 = calc.find_homologous(z, &r, false).unwrap().expect("expected Some for trefoil R(0,1)");
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

    // The expected regions change only the order: a right hint and a wrong one give the same G0.
    #[test]
    fn test_trefoil_expected() {
        let gd = GridDiagram::test_data("3_1").unwrap();
        for hint in [vec![Region::simple(0, 1), Region::simple(1, 0)], vec![Region::simple(0, 0)]] {
            let mut calc = G0Calc::init(&gd, None, None, None, config());
            calc.set_expected(hint);
            calc.run();

            let rs = calc.accepted();
            assert_eq!(rs.len(), 2);
            assert!(rs.contains(&Region::simple(0, 1)));
            assert!(rs.contains(&Region::simple(1, 0)));
        }
    }

    // A resigned expected region leaves the candidates (but not `expected`), so the run ends.
    // note: R(0, 1) and R(1, 0) contain z, so they need no system. R(0, 0) needs one, which
    // is tiny for the trefoil, so max_size = 0 to make it resign.
    #[test]
    fn test_resigned_expected_terminates() {
        let gd = GridDiagram::test_data("3_1").unwrap();
        let cfg = G0CalcConfig { max_size: 0, ..config() };
        let mut calc = G0Calc::init(&gd, None, None, None, cfg);
        calc.set_expected(vec![Region::simple(0, 0)]);
        calc.run();

        assert!(calc.candidates.regions.is_empty());
        assert!(calc.resigned().contains(&Region::simple(0, 0)),
            "accepted {:?}, resigned {:?}", calc.accepted(), calc.resigned());
        assert!(calc.accepted().contains(&Region::simple(0, 1)));
        assert_eq!(calc.candidates.expected, vec![Region::simple(0, 0)]);
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
