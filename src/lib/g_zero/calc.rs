use log::*;
use num_traits::{One, Zero};
use yui_core::lc::{Gen, Lc};
use yui_core::num::FF2;
use yui_core::poly::MultiVar;
use yui_homology::{DisplayTable, Summand};
use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};
use yui_matrix::sparse::pluq::{solve_pluq, solve_pluq_incr};

use crate::{GCGen, GridComplex, GridDiagram, HFK};
use crate::grid_complex::{GCFiltered, GCInfGen};
use crate::grid_complex::chain::GCInfChainExt;
use super::{Candidates, Region};

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct G0CalcConfig {
    pub filter_regions: bool,
    pub filter_gens: bool,
    pub accept_if_unique: bool,
    pub try_fast_ng: bool,
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
            max_size: 1_000_000,
            max_piv: 500_000,
            chunk: 100_000,
        }
    }
}

/// Returned by `find_homologous` when the system exceeds `max_size` and
/// neither the trivial cycle-membership check nor the U-specialization
/// fast-path concluded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum G0CalcError { Resigned }

pub struct G0Calc {
    gc: GridComplex<GCFiltered>,
    g: usize,
    z: Lc<GCInfGen, FF2>,
    candidates: Candidates,
    config: G0CalcConfig
}

impl G0Calc {
    pub(crate) fn new(
        gc: GridComplex<GCFiltered>,
        g: usize,
        z: Lc<GCInfGen, FF2>,
        candidates: Candidates,
        config: G0CalcConfig,
    ) -> Self {
        G0Calc { gc, g, z, candidates, config }
    }

    pub fn from_gd(gd: &GridDiagram, g4: Option<usize>, tau: Option<isize>, config: G0CalcConfig) -> Self {
        info!("compute HFK ...");
        let hfk = HFK::from(gd);
        info!("HFK:\n{}", hfk.display_table());

        let g = hfk.genus() as usize;
        info!("g = {g}");

        Self::validate(&hfk, g, g4, tau);

        let candidates = Self::make_candidates(&hfk, g4, tau, config);

        info!("prepare generators ...");
        let gc = Self::build_complex(gd, &candidates, config.filter_gens);
        info!("  #generators: {}", gc.gens.len());

        info!("compute homological generator..");
        let z = gc.hom_gen();
        
        if z.is_zero() { 
            info!("  z = 0.");
        } else {
            info!("  support of z: \n{}", z.support_grid().display_table());
        }

        Self::new(gc, g, z, candidates, config)
    }

    fn validate(hfk: &HFK, g: usize, g4: Option<usize>, tau: Option<isize>) { 
        // MEMO: replace with if-let chain
        if let Some(g4) = g4 { if let Some(tau) = tau { 
            let t_abs = tau.abs() as usize;
            assert!(t_abs <= g4, "Invalid inequality: |tau| = {t_abs} > g4 = {g4}");
        }}
        
        if let Some(g4) = g4 { 
            assert!(g4 <= g, "Invalid inequality: g4 = {g4} > g = {g}");
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
            candidates.deduplicate();
            debug!("  deduplicated -> {}", candidates.regions.len());

            candidates.filter_by_hfk(&hfk);
            debug!("  filtered by hfk -> {}", candidates.regions.len());

            if let Some(g4) = g4 { 
                candidates.filter_by_g4(g4);
                debug!("  filtered by g4({g4}) -> {}", candidates.regions.len());
            }

            if let Some(tau) = tau { 
                candidates.filter_by_tau(tau);
                debug!("  filtered by tau({tau}) -> {}", candidates.regions.len());
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

    pub(crate) fn build_complex(gd: &GridDiagram, candidates: &Candidates, filter_gens: bool) -> GridComplex<GCFiltered> {
        let gd = gd.clone();
        if !filter_gens { return GridComplex::new(gd); }

        GridComplex::new_filtered(gd, |x| 
            candidates.regions.iter().any(|r| { 
                let s = r.right_end(); // r is to be shifted by (-s, -s).
                let m = x.maslov_deg;
                if m < -2 * s { return false; }

                // Consider U^e x such that 
                //
                //    M(U^e x) = -2s or -2s + 1.
                //
                // Then, 
                //
                //   A(U^e x)   = -e + A(x). 
                //   alg(U^e x) = -e. 
                //
                // Check U^e x ∈ R[-s, -s]

                let e = (m + 2 * s) / 2; // >= 0
                !r.contains((s - e, s - e + x.alex_deg)) // <=> !r.shift(-s, -s).contains(-e, -e + x.alex_deg)
            })
        )
    }

    pub fn complex(&self) -> &GridComplex<GCFiltered> {
        &self.gc
    }

    pub fn genus(&self) -> usize {
        self.g
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
            self.candidates.accept(&r);
            return
        }

        let s = r.right_end().max(0);
        let r_shifted = r.shift(-s, -s);
        let z_shifted = self.shift_z(s as usize);

        match self.find_homologous(&z_shifted, &r_shifted) {
            Ok(Some(_)) => {
                info!("OK: {r}");
                self.candidates.accept(&r);
            }
            Ok(None) => {
                info!("NG: {r}");
                self.candidates.reject(&r);
            }
            Err(G0CalcError::Resigned) => {
                info!("RESIGN: {r}");
                self.candidates.resign(&r);
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

        // try NG fast-path for large system
        if self.config.try_fast_ng && (n_w >= self.config.chunk || n_z >= self.config.chunk) {
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

        let Some(x_sol) = self.find_sol(&a, &y) else {
            return Ok(None);
        };
        let w = gen_w.devectorize(&x_sol);

        // z - d(w) is homologous to z and lies inside r.
        let mut z1 = z.clone();
        z1 -= &self.gc.differentiate(&w);
        Ok(Some(z1))
    }

    fn find_sol(&self, a: &SpMat<FF2>, y: &SpVec<FF2>) -> Option<SpVec<FF2>> {
        debug!("solve linear system ...");

        let x = if a.nrows() < self.config.chunk {
            solve_pluq(a, y)
        } else {
            solve_pluq_incr(a, y, self.config.max_piv, self.config.chunk)
        }?;

        debug!("validate solution..");
        assert_eq!(&(a * &x), y);
        debug!("  ok.");

        Some(x)
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

        if !self.has_sol_specialed(&raw_gen_w, &raw_gen_z, z0, r, &vec![FF2::one(); n]) {
            debug!("  NG via U = 1");
            return false;
        }

        debug!("passed NG fast-path.");
        true
    }

    fn has_sol_specialed(
        &self,
        raw_gen_w: &Summand<GCGen, FF2>,
        raw_gen_z: &Summand<GCGen, FF2>,
        z0: &Lc<GCInfGen, FF2>,
        r: &Region,
        c: &[FF2],
    ) -> bool {
        let lift = |x: &GCGen| Lc::from(GCInfGen::from(x.clone()));
        let proj = move |x: &GCInfGen| if x.u.eval::<FF2>(c).is_zero() { 
            Lc::zero() 
        } else { 
            Lc::from(x.gen.clone()) 
        };

        debug!("create specialized matrix ...");
        let (a_c, y_c) = create_ay(
            raw_gen_w, raw_gen_z, &z0.apply(&proj),
            |x| r.cutoff(&self.gc.differentiate(&x.apply(&lift))).apply(&proj),
        );

        debug!("  specialized: {:?} (nnz = {})", a_c.shape(), a_c.nnz());
        self.find_sol(&a_c, &y_c).is_some()
    }
}

fn create_ay<X>(
    gen_w: &Summand<X, FF2>,
    gen_z: &Summand<X, FF2>,
    y_target: &Lc<X, FF2>,
    f: impl Fn(&Lc<X, FF2>) -> Lc<X, FF2>,
) -> (SpMat<FF2>, SpVec<FF2>)
where
    X: Gen + Into<GCInfGen>,
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
        let calc = G0Calc::from_gd(&gd, None, None, cfg);

        let z = &calc.z;
        assert!(!z.is_zero());

        let r = Region::simple(0, 0);
        assert!(r.contains_chain(&z));
        
        let z1 = calc.find_homologous(&z, &r).unwrap().expect("expected Some");
        assert!(z1.is_zero());
    }

    #[test]
    fn test_find_homologous_unknot_none() {
        // R[(-1,0), (0,-1)] lies strictly below the hom_gen at (0,0).
        let gd = GridDiagram::unknot();
        let cfg = config();
        let calc = G0Calc::from_gd(&gd, None, None, cfg);
        let z = &calc.z;

        assert!(!z.is_zero());

        let r = Region::new([Corner::new(-1, 0), Corner::new(0, -1)]);
        assert_eq!(calc.find_homologous(&z, &r), Ok(None));
    }

    #[test]
    fn test_find_homologous_trefoil() {
        // For the trefoil (genus 1), R(0,1) is the expected realizable region.
        // find_homologous should return Some(z1) where z1 is a cycle inside R(0,1).
        let gd = GridDiagram::trefoil();
        let cfg = config();
        let calc = G0Calc::from_gd(&gd, None, None, cfg);
        let z = &calc.z;
        assert!(!z.is_zero());

        let r = Region::simple(0, 1);
        assert!(r.contains_chain(&z));

        let z1 = calc.find_homologous(&z, &r).unwrap().expect("expected Some for trefoil R(0,1)");
        assert!(z1.is_zero());
    }

    fn accepted(gd: GridDiagram) -> Vec<Region> {
        let cfg = config();
        let mut calc = G0Calc::from_gd(&gd, None, None, cfg);
        calc.run();

        let mut rs = calc.accepted().to_vec();
        rs.sort_by_key(|r| r.right_end());
        rs
    }

    #[test]
    fn test_unknot() {
        let rs = accepted(GridDiagram::load("0_1").unwrap());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_trefoil() {
        let rs = accepted(GridDiagram::load("3_1").unwrap());
        assert_eq!(rs.len(), 2);
        assert!(rs.contains(&Region::simple(0, 1)));
        assert!(rs.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_trefoil_mirror() {
        let rs = accepted(GridDiagram::load("3_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([Corner::new(-1, 0), Corner::new(0, -1)])]);
    }

    #[test]
    fn test_figure8() {
        let rs = accepted(GridDiagram::load("4_1").unwrap());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_figure8_mirror() {
        let rs = accepted(GridDiagram::load("4_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_5_1() {
        let rs = accepted(GridDiagram::load("5_1").unwrap());
        assert_eq!(rs.len(), 3);
        assert!(rs.contains(&Region::simple(0, 2)));
        assert!(rs.contains(&Region::simple(1, 1)));
        assert!(rs.contains(&Region::simple(2, 0)));
    }

    #[test]
    fn test_5_1_mirror() {
        let rs = accepted(GridDiagram::load("5_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([
            Corner::new(-2, 0),
            Corner::new(-1, -1),
            Corner::new(0, -2),
        ])]);
    }

    #[test]
    fn test_5_2() {
        let rs = accepted(GridDiagram::load("5_2").unwrap());
        assert_eq!(rs.len(), 2);
        assert!(rs.contains(&Region::simple(0, 1)));
        assert!(rs.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_5_2_mirror() {
        let rs = accepted(GridDiagram::load("5_2").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([Corner::new(-1, 0), Corner::new(0, -1)])]);
    }
}
