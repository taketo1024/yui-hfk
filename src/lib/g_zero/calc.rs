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

pub struct G0CalcConfig {
    pub filter_regions: bool,
    pub filter_gens: bool,
    pub max_piv: usize,
    pub chunk: usize,
}

impl Default for G0CalcConfig {
    fn default() -> Self {
        Self {
            filter_regions: true,
            filter_gens: true,
            max_piv: 500_000,
            chunk: 100_000,
        }
    }
}

pub struct G0Calc {
    gc: GridComplex<GCFiltered>,
    z: Lc<GCInfGen, FF2>,
    candidates: Candidates,
    max_piv: usize,
    chunk: usize,
}

impl G0Calc {
    pub fn new(gd: GridDiagram, config: G0CalcConfig) -> Self {
        info!("compute HFK ...");
        let hfk = HFK::from(&gd);
        info!("HFK:\n{}", hfk.display_table());

        let g = hfk.genus() as usize;
        info!("g = {g}");

        let candidates = Self::make_candidates(&hfk, g, config.filter_regions);
        info!("candidate regions: {}", candidates.regions.len());
        for r in &candidates.regions {
            info!("  {r}");
        }

        info!("prepare generators ...");
        let gc = Self::build_complex(gd, &candidates, config.filter_gens);
        info!("  #generators: {}", gc.gens.len());

        info!("compute homological generator");
        let z = gc.hom_gen();

        if z.is_zero() { 
            info!("  homological generator contained in all candidate regions!");
        }

        G0Calc { gc, z, candidates, max_piv: config.max_piv, chunk: config.chunk }
    }

    fn make_candidates(hfk: &HFK, g: usize, filter_regions: bool) -> Candidates { 
        let mut candidates = Candidates::generate_from(g);
        if filter_regions { 
            candidates.deduplicate();
            candidates.filter_by_hfk(&hfk);
        }
        candidates.sort();
        candidates
    }

    fn build_complex(gd: GridDiagram, candidates: &Candidates, filter_gens: bool) -> GridComplex<GCFiltered> {
        if !filter_gens { return GridComplex::new(gd); }

        let s_max = candidates.regions.iter()
            .map(|r| r.right_end().max(0))
            .max()
            .unwrap_or(0);
        let s_range = -2 * s_max ..= 1;

        GridComplex::new_filtered(gd, |x| {
            if !s_range.contains(&x.maslov_deg) { return false; }
            let s = -(x.maslov_deg.div_euclid(2)); // >= 0
            candidates.regions.iter()
                .filter(|r| r.right_end() == s)
                .any(|r| !r.contains((s, x.alex_deg + s))) // equivalent to: !r.shift(-s, -s).contains(0, x.alex_deg)
        })
    }

    pub fn complex(&self) -> &GridComplex<GCFiltered> { 
        &self.gc
    }

    /// Processes all remaining candidates and populates `candidates.accepted`.
    pub fn run(&mut self) {
        while !self.candidates.regions.is_empty() {
            let r = self.candidates.regions.remove(0);
            self.process(r);
        }
    }

    pub fn filter_by_g4(&mut self, g4: usize) {
        self.candidates.filter_by_g4(g4);
        info!("g4 = {g4}, regions: {}", self.candidates.regions.len());
    }

    pub fn accepted(&self) -> &[Region] {
        &self.candidates.accepted
    }

    // Tests whether `r` is realizable: tries to find a cycle homologous to z
    // inside r. Shifts both r and z by s = right_end(r) so that the shifted
    // region has right_end = 0, which is required by find_homologous.
    fn process(&mut self, r: Region) {
        if r.contains_chain(&self.z) {
            self.candidates.accept(&r);
            return
        }

        info!("try: {r} (remain: {})", self.candidates.regions.len());

        let s = r.right_end().max(0);
        let r_shifted = r.shift(-s, -s);
        let z_shifted = self.shift_z(s as usize);
        match self.find_homologous(&z_shifted, &r_shifted) {
            Some(_) => { 
                info!("OK: {r}"); 
                self.candidates.accept(&r); 
            },
            None => { 
                info!("NG: {r}"); 
                self.candidates.reject(&r); 
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
    fn find_homologous(
        &self,
        z: &Lc<GCInfGen, FF2>,
        r: &Region,
    ) -> Option<Lc<GCInfGen, FF2>> {
        assert!(z.is_maslov_homogeneous());

        let k = z.maslov_degree();
        let z0 = r.cutoff(z);

        let pts_w = self.gc.complementary_points(k + 1, r);
        let pts_z = self.gc.complementary_points(k, r);

        let n_w = self.gc.count_inf_gens_in(k + 1, &pts_w);
        let n_z = self.gc.count_inf_gens_in(k, &pts_z);

        debug!("linear system size: {:?}", (n_z, n_w));

        // try NG fast-path for large system
        if n_w >= self.chunk || n_z >= self.chunk {
            let n_raw_w = self.gc.count_raw_gens_in(k + 1, &pts_w);
            let n_raw_z = self.gc.count_raw_gens_in(k, &pts_z);
            if (n_raw_w < n_w || n_raw_z < n_z) && !self.pass_specialized(k, &pts_w, &pts_z, &z0, r) {
                return None;
            }
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

        let x_sol = self.find_sol(&a, &y)?;
        let w = gen_w.devectorize(&x_sol);

        // z - d(w) is homologous to z and lies inside r.
        let mut z1 = z.clone();
        z1 -= &self.gc.differentiate(&w);
        Some(z1)
    }

    fn find_sol(&self, a: &SpMat<FF2>, y: &SpVec<FF2>) -> Option<SpVec<FF2>> {
        debug!("solve linear system ...");

        let y = y.to_dense();

        let x = if a.nrows() < self.chunk {
            solve_pluq(a, &y)
        } else {
            solve_pluq_incr(&a, &y, self.max_piv, self.chunk)
        }?;

        Some(SpVec::from(x))
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

        if !self.has_sol_specialed(&raw_gen_w, &raw_gen_z, z0, r, &vec![FF2::zero(); n]) {
            debug!("  NG via U = 0");
            return false;
        }

        if !self.has_sol_specialed(&raw_gen_w, &raw_gen_z, z0, r, &vec![FF2::one(); n]) {
            debug!("  NG via U = 1");
            return false;
        }

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
        let proj = move |x: &GCInfGen| if x.u.eval(c).is_zero() { 
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

    #[test]
    fn test_find_homologous_unknot_already_inside() {
        // For the unknot the hom_gen lies at (alg=0, alex=0) ∈ R(0,0),
        // so z0 = cutoff(z, R(0,0)) = 0.  The system is trivially satisfied
        // with w = 0, and the function returns z unchanged.
        let gd = GridDiagram::unknot();
        let cfg = G0CalcConfig { filter_regions: false, filter_gens: false, ..Default::default() };
        let calc = G0Calc::new(gd, cfg);
        let cx = calc.complex();
        let z = cx.hom_gen();

        assert!(!z.is_zero());

        let r = Region::simple(0, 0);
        let z1 = calc.find_homologous(&z, &r).expect("expected Some");

        assert!(!z1.is_zero());
        assert!(cx.differentiate(&z1).is_zero());
        assert!(r.contains_chain(&z1));
    }

    #[test]
    fn test_find_homologous_unknot_none() {
        // R[(-1,0), (0,-1)] lies strictly below the hom_gen at (0,0).
        let gd = GridDiagram::unknot();
        let cfg = G0CalcConfig { filter_regions: false, filter_gens: false, ..Default::default() };
        let calc = G0Calc::new(gd, cfg);
        let cx = calc.complex();
        let z = cx.hom_gen();
        
        assert!(!z.is_zero());

        let r = Region::new([Corner::new(-1, 0), Corner::new(0, -1)]);
        let z1 = calc.find_homologous(&z, &r);

        assert!(z1.is_none());
    }

    #[test]
    fn test_find_homologous_trefoil() {
        // For the trefoil (genus 1), R(0,1) is the expected realizable region.
        // find_homologous should return Some(z1) where z1 is a cycle inside R(0,1).
        let gd = GridDiagram::trefoil();
        let cfg = G0CalcConfig { filter_regions: false, filter_gens: false, ..Default::default() };
        let calc = G0Calc::new(gd, cfg);
        let cx = calc.complex();
        let z = cx.hom_gen();
        assert!(!z.is_zero());

        let r = Region::simple(0, 1);
        let z1 = calc.find_homologous(&z, &r).expect("expected Some for trefoil R(0,1)");

        assert!(!z1.is_zero());
        assert!(cx.differentiate(&z1).is_zero(), "z1 must be a cycle");
        assert!(r.contains_chain(&z1), "z1 must lie inside r");
    }

    fn accepted(gd: GridDiagram) -> Vec<Region> {
        let cfg = G0CalcConfig { filter_regions: false, filter_gens: false, ..Default::default() };
        let mut calc = G0Calc::new(gd, cfg);
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
