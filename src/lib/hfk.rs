use std::collections::HashMap;
use std::ops::Index;

use delegate::delegate;
use num_integer::binomial;
use num_traits::Zero;
use yui_homology::{DisplayTable, Grid2, GridTrait, SummandTrait, isize2, range_of};

use crate::{GCTilde, GridComplex, GridDiagram};

/// The bigraded module structure of ĤFK(K).
pub struct HFK {
    data: Grid2<usize>,
}

impl HFK {
    /// The 3-genus of K, equal to the maximum Alexander degree.
    pub fn genus(&self) -> isize {
        self.data.support().map(|i| i.1).max().unwrap_or(0)
    }

    /// Mirror image: negate both Maslov and Alexander degrees.
    pub fn mirror(&self) -> Self {
        let support: Vec<isize2> = self.data.support().map(|i| isize2(-i.0, -i.1)).collect();
        let data = Grid2::generate(support, |i| *self.data.get(isize2(-i.0, -i.1)));
        HFK { data }
    }

    delegate! { 
        to self.data { 
            pub fn iter(&self) -> impl Iterator<Item = (isize2, &usize)>;
            pub fn display_table(&self, l1: &str, l2: &str) -> String;
            pub fn print_table(&self, l1: &str, l2: &str);
        }
    }
}

impl From<&GridDiagram> for HFK {
    fn from(gd: &GridDiagram) -> Self {
        let n = gd.grid_number();

        let cx = GridComplex::<GCTilde>::new_filtered(gd.clone(), |x| x.alex_deg >= 0);
        let m_range = cx.maslov_deg_range();
        let a_range = cx.alex_deg_range();
        let h = cx.into_homology();

        let mut ranks: HashMap<(isize, isize), usize> = HashMap::new();

        for j in a_range.clone() { 
            for i in m_range.clone() { 
                if !h[(i, j)].is_zero() { 
                    ranks.insert((i, j), h[(i, j)].rank());
                }
            }
        }

        // divide by W^(n - 1).
        for j in a_range.clone().rev() { 
            for i in m_range.clone() { 
                let Some(&r) = ranks.get(&(i, j)) else { continue };
                if r.is_zero() { continue }

                for k in 1 .. n { 
                    let k = k as isize;
                    if k > j { break }
                    
                    let t = ranks.get_mut(&(i - k, j - k)).unwrap();
                    *t -= r * binomial(n - 1, k as u8) as usize;
                }
            }
        }

        ranks.retain(|_, r| *r > 0);

        // use symmetry: HFK_{i, j}(K) = HFK_{i - 2j, -j}(K)
        for j in a_range.clone().rev() { 
            if j == 0 { break }

            for i in m_range.clone() { 
                let Some(r) = ranks.get(&(i, j)) else { continue };
                ranks.insert((i - 2 * j, -j), *r);
            }
        }

        let m_range = range_of(ranks.keys().map(|&(i, _)| i));
        let a_range = range_of(ranks.keys().map(|&(_, j)| j));
        let support = m_range.flat_map(|i| a_range.clone().map(move |j| isize2(i, j)));

        let data = Grid2::generate(
            support,
            |isize2(i, j)| *ranks.get(&(i, j)).unwrap_or(&0)
        );

        HFK { data }
    }
}

impl FromIterator<(isize2, usize)> for HFK {
    fn from_iter<T: IntoIterator<Item = (isize2, usize)>>(iter: T) -> Self {
        let data = Grid2::from_iter(iter);
        Self { data }
    }
}

impl Index<(isize, isize)> for HFK {
    type Output = usize;
    fn index(&self, index: (isize, isize)) -> &Self::Output {
        &self.data[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hfk_unknot() {
        let gd = GridDiagram::unknot();
        let h = HFK::from(&gd);
        assert_eq!(h[(0, 0)], 1);
    }

    #[test]
    fn test_hfk_trefoil() {
        let gd = GridDiagram::trefoil();
        let h = HFK::from(&gd);
        assert_eq!(h[(0, 1)], 1);
        assert_eq!(h[(-1, 0)], 1);
        assert_eq!(h[(-2, -1)], 1);
    }

    #[test]
    fn test_hfk_trefoil_mirror() {
        let gd = GridDiagram::trefoil().rotate90();
        let h = HFK::from(&gd);
        assert_eq!(h[(2, 1)], 1);
        assert_eq!(h[(1, 0)], 1);
        assert_eq!(h[(0, -1)], 1);
    }

    #[test]
    fn test_hfk_4_1() {
        let gd = GridDiagram::figure8();
        let h = HFK::from(&gd);
        assert_eq!(h[(1, 1)], 1);
        assert_eq!(h[(0, 0)], 3);
        assert_eq!(h[(-1, -1)], 1);
    }
}