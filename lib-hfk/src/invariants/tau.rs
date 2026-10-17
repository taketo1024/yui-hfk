use log::{debug, info};

use crate::g_zero::{Candidates, Corner, Region};
use crate::g_zero::calc::G0CalcError;
use crate::{G0Calc, G0CalcConfig, GridDiagram, genus};

pub fn tau(gd: &GridDiagram) -> Result<isize, G0CalcError> { 
    tau_with_config(gd, None, None, G0CalcConfig::default())
}

pub fn tau_with_config(gd: &GridDiagram, g: Option<usize>, g4: Option<usize>, config: G0CalcConfig) -> Result<isize, G0CalcError> {
    info!("compute tau ...");

    let g = g.unwrap_or_else(|| { 
        debug!("compute genus ...");
        genus(gd)
    });
    debug!("  g = {g}");

    if let Some(g4) = g4 { 
        assert!(g4 <= g, "Invalid inequality: g4 = {g4} > g = {g}");
    }

    let t = match compute_tau(gd, g, g4, config) {
        Ok(t) if t >= 0 => t,
        Ok(_) => {
            info!("reached tau <= -1. switch to mirror.");
            let t = compute_tau(&gd.rotate90(), g, g4, config)?;
            assert!(t >= 1);
            -t
        }
        Err(G0CalcError::Resigned) => {
            info!("resigned. switch to mirror.");
            match compute_tau(&gd.rotate90(), g, g4, config)? {
                t if t >= 0 => -t,
                // note: tau(mirror) ≤ -1 only gives tau ≥ 1.
                _ => return Err(G0CalcError::Resigned),
            }
        }
    };

    info!("computed tau = {t}");
    validate(g4, t);
    Ok(t)
}

/// `τ(K) = min{ m ∈ ℤ | ∃R ∈ G0(K), R ⊂ {i ≤ −1} ∪ R(0, m) }` ([SS] Thm 2.20), where `R` is
/// included iff every corner `(i, j)` of `R` has `i ≤ −1`, or `i = 0` and `j ≤ m`.
pub fn tau_from_g0(g0: &[Region]) -> isize {
    g0.iter()
        .filter(|r| r.corners().iter().all(|c| c.i <= 0))
        .map(|r| r.corners().iter().filter(|c| c.i == 0).map(|c| c.j).max().unwrap_or(isize::MIN))
        .min()
        .expect("a region of G0 lies in {i ≤ 0}")
}

fn compute_tau(gd: &GridDiagram, g: usize, g4: Option<usize>, config: G0CalcConfig) -> Result<isize, G0CalcError> {
    let mut candidates = Candidates::new();
    let t = g4.unwrap_or(g) as isize;

    // MEMO: Only check regions r for i >= -1. 
    for i in (-1..=t).rev() { 
        let r = Region::new([Corner::new(-1, isize::MAX), Corner::new(0, i)]);
        candidates.add(r);
    }

    let mut calc = G0Calc::init_with_candidates(gd, candidates, config);
    calc.run();

    if !calc.resigned().is_empty() { 
        return Err(G0CalcError::Resigned)
    }

    if calc.accepted().len() == 2 { 
        let r = calc.accepted().iter().find(|r| r.corners()[1].i == 0).unwrap();
        let t = r.corners()[1].j;
        Ok(t)
    } else { 
        panic!("Invalid state, possibly due to invalid input: g = {g}, g4 = {g4:?}");
    }
}

fn validate(g4: Option<usize>, tau: isize) { 
    if let Some(g4) = g4 {
        let t_abs = tau.unsigned_abs();
        assert!(t_abs <= g4, "Invalid inequality: |tau| = {t_abs} > g4 = {g4}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::g_zero::dual;

    #[test]
    fn test_unknot() {
        assert_eq!(tau(&GridDiagram::unknot()), Ok(0));
    }

    #[test]
    fn test_trefoil() {
        assert_eq!(tau(&GridDiagram::trefoil()), Ok(1));
    }

    #[test]
    fn test_figure8() {
        assert_eq!(tau(&GridDiagram::figure8()), Ok(0));
    }

    #[test]
    fn test_5_1() {
        assert_eq!(tau(&GridDiagram::test_data("5_1").unwrap()), Ok(2));
    }

    #[test]
    fn test_5_2() {
        assert_eq!(tau(&GridDiagram::test_data("5_2").unwrap()), Ok(1));
    }

    #[test]
    fn test_6_1() {
        assert_eq!(tau(&GridDiagram::test_data("6_1").unwrap()), Ok(0));
    }

    #[test]
    fn test_6_2() {
        assert_eq!(tau(&GridDiagram::test_data("6_2").unwrap()), Ok(1));
    }

    #[test]
    fn test_6_3() {
        assert_eq!(tau(&GridDiagram::test_data("6_3").unwrap()), Ok(0));
    }

    // tau is antisymmetric under mirroring.
    #[test]
    fn test_mirror_antisymmetry() {
        for name in ["3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3"] {
            let gd = GridDiagram::test_data(name).unwrap();
            let t = tau(&gd).unwrap();
            let t_mirror = tau(&gd.rotate90()).unwrap();
            assert_eq!(t, -t_mirror, "{name}: tau = {t}, mirror = {t_mirror}");
        }
    }

    #[test]
    fn test_unknot_mirror() {
        assert_eq!(tau(&GridDiagram::unknot().rotate90()), Ok(0));
    }

    #[test]
    fn test_trefoil_mirror() {
        assert_eq!(tau(&GridDiagram::trefoil().rotate90()), Ok(-1));
    }

    #[test]
    fn test_figure8_mirror() {
        assert_eq!(tau(&GridDiagram::figure8().rotate90()), Ok(0));
    }

    // note: for the mirror of 6_1, the largest system has size 247, and 62 on its mirror.
    #[test]
    fn test_resigned_switches_to_mirror() {
        let gd = GridDiagram::test_data("6_1").unwrap().rotate90();
        let config = |max_size| G0CalcConfig { max_size, ..G0CalcConfig::default() };
        assert_eq!(compute_tau(&gd, genus(&gd), None, config(100)), Err(G0CalcError::Resigned));
        assert_eq!(tau_with_config(&gd, None, None, config(100)), Ok(0));
        assert_eq!(tau_with_config(&gd, None, None, config(20)), Err(G0CalcError::Resigned));
    }

    fn g0(ty: &[&[(isize, isize)]]) -> Vec<Region> {
        ty.iter().map(|r| Region::new(r.iter().map(|&(i, j)| Corner::new(i, j)))).collect()
    }

    // τ(K) of Table 2 in [SS] (v2), and τ(K*) = −τ(K) with G0(K*) = dual(G0(K)).
    #[test]
    fn test_tau_from_g0_table_2() {
        let table: [(&[&[(isize, isize)]], isize); 7] = [
            (&[&[(0, 0)]], 0),
            (&[&[(0, 1)], &[(1, 0)]], 1),
            (&[&[(0, 2)], &[(1, 1)], &[(2, 0)]], 2),
            (&[&[(0, 2)], &[(0, 1), (1, 0)], &[(2, 0)]], 2),
            (&[&[(0, 3)], &[(1, 2)], &[(2, 1)], &[(3, 0)]], 3),
            (&[&[(0, 3)], &[(1, 1)], &[(3, 0)]], 3),
            (&[&[(0, 4)], &[(1, 2)], &[(2, 1)], &[(4, 0)]], 4),
        ];
        for (ty, t) in table {
            assert_eq!(tau_from_g0(&g0(ty)), t, "{ty:?}");
            assert_eq!(tau_from_g0(&dual(&g0(ty))), -t, "{ty:?}: mirror");
        }
    }

    #[test]
    fn test_tau_from_g0_5_1() {
        let gd = GridDiagram::test_data("5_1").unwrap();
        let mut calc = G0Calc::init(&gd, None, None, None, G0CalcConfig::default());
        calc.run();
        assert_eq!(Ok(tau_from_g0(&calc.g0().unwrap())), tau(&gd));
    }
}