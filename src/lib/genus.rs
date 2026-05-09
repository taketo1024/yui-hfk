use yui_homology::SummandTrait;

use crate::{GCTilde, GridComplex, GridDiagram};

pub fn genus(gd: &GridDiagram) -> usize { 
    let cx = GridComplex::<GCTilde>::new_filtered(gd.clone(), |x| x.alex_deg > 0);
    let m_range = cx.maslov_deg_range();
    let a_range = cx.alex_deg_range();

    let h = cx.into_homology();
    a_range.rev()
        .find_map(|j| {
            m_range.clone()
                .any(|i| h[(i, j)].rank() > 0)
                .then_some(j as usize)
        })
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genus_unknot() {
        let gd = GridDiagram::unknot();
        let g = genus(&gd);
        assert_eq!(g, 0);
    }

    #[test]
    fn test_genus_trefoil() {
        let gd = GridDiagram::trefoil();
        let g = genus(&gd);
        assert_eq!(g, 1);
    }

    #[test]
    fn test_genus_5_1() {
        let gd = GridDiagram::from_sequences(&[5,4,3,2,1,0,6], &[0,6,5,4,3,2,1]);
        let g = genus(&gd);
        assert_eq!(g, 2);
    }
}