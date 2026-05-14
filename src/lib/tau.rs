use log::{debug, info};

use crate::g_zero::{Candidates, Corner, Region};
use crate::g_zero::calc::Resigned;
use crate::{G0Calc, G0CalcConfig, GridDiagram, genus};

pub fn tau(gd: &GridDiagram) -> Result<isize, Resigned> { 
    tau_with_config(gd, None, None, G0CalcConfig::default())
}

pub fn tau_with_config(gd: &GridDiagram, g: Option<usize>, g4: Option<usize>, config: G0CalcConfig) -> Result<isize, Resigned> { 
    info!("compute tau ...");

    let g = g.unwrap_or_else(|| { 
        debug!("compute genus ...");
        genus(&gd)
    });
    debug!("  g = {g}");

    let t = compute_tau(gd, g, g4, config)?;

    if t > 0 { 
        Ok(t) // minimal 
    } else { 
        let t = compute_tau(&gd.rotate90(), g, g4, config)?;
        Ok(-t)
    }
}

fn compute_tau(gd: &GridDiagram, g: usize, g4: Option<usize>, config: G0CalcConfig) -> Result<isize, Resigned> { 
    let mut candidates = Candidates::new();
    let t = g4.unwrap_or(g) as isize;

    // MEMO: Only check regions r for i >= 0. 
    for i in (0..=t).rev() { 
        let r = Region::new([Corner::new(-1, isize::MAX), Corner::new(0, i)]);
        candidates.add(r);
    }

    info!("prepare generators ...");
    let gc = G0Calc::build_complex(gd, &candidates, true);
    info!("  #generators: {}", gc.gens.len());

    info!("compute homological generator");
    let z = gc.hom_gen();
    let g = g as usize;
    let mut calc = G0Calc::new(gc, g, z, candidates, config);

    calc.run();

    if !calc.resigned().is_empty() { 
        return Err(Resigned)
    }

    assert_eq!(calc.accepted().len(), 2);

    let r = calc.accepted().iter().find(|r| r.corners()[1].i == 0).unwrap();
    let t = r.corners()[1].j;
    Ok(t)
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