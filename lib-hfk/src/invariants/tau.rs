use log::{debug, info};

use crate::g_zero::{Candidates, Corner, Region};
use crate::g_zero::calc::G0CalcError;
use crate::{G0Calc, G0CalcConfig, GridDiagram, genus};

pub fn tau(gd: &GridDiagram) -> Result<isize, G0CalcError> { 
    tau_with_config(gd, None, None, G0CalcConfig::default())
}

pub fn tau_with_config(gd: &GridDiagram, g: Option<usize>, g4: Option<usize>, config: G0CalcConfig) -> Result<isize, G0CalcError> {
    // The search is markedly cheaper on the positive side, and tau(mirror) =
    // -tau, so start from there. Mirroring negates the writhe, so this recurses
    // at most once.
    if gd.writhe() < 0 {
        info!("negative writhe, switch to mirror.");
        return tau_with_config(&gd.rotate90(), g, g4, config).map(|t| -t);
    }

    info!("compute tau ...");

    let g = g.unwrap_or_else(|| { 
        debug!("compute genus ...");
        genus(gd)
    });
    debug!("  g = {g}");

    if let Some(g4) = g4 { 
        assert!(g4 <= g, "Invalid inequality: g4 = {g4} > g = {g}");
    }

    let t = compute_tau(gd, g, g4, config)?;

    if t >= 0 { 
        info!("computed tau = {}", t);
        validate(g4, t);
        Ok(t) // minimal 
    } else if t == -1 { 
        info!("reached tau <= -1. switch to mirror.");
        let t = compute_tau(&gd.rotate90(), g, g4, config)?;

        assert!(t >= 1);
        info!("computed tau = {}", -t);

        validate(g4, -t);
        Ok(-t)
    } else { 
        panic!("Unreachable.")
    }
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
}