use itertools::Itertools;

use super::{Corner, Region};

/// `G0(K*)` from `G0(K)` by [SS] Cor. A.7: the minimal regions among `cl({-p₁, …, -p_N})`,
/// over all choices of a corner `pᵢ` of each region `Rᵢ ∈ G0(K)`.
pub fn dual(g0: &[Region]) -> Vec<Region> {
    let candidates: Vec<Region> = g0.iter()
        .map(|r| r.corners().iter().copied())
        .multi_cartesian_product()
        .map(|ps| Region::new(ps.into_iter().map(|p| Corner::new(-p.i, -p.j))))
        .unique()
        .collect();

    candidates.iter()
        .filter(|r| !candidates.iter().any(|s| s != *r && r.includes(s)))
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    fn regions(ty: &[&[(isize, isize)]]) -> HashSet<Region> {
        ty.iter().map(|r| Region::new(r.iter().map(|&(i, j)| Corner::new(i, j)))).collect()
    }

    fn dual_set(g0: &HashSet<Region>) -> HashSet<Region> {
        dual(&g0.iter().cloned().collect_vec()).into_iter().collect()
    }

    #[test]
    fn dual_of_t2_3() {
        let g0 = regions(&[&[(0, 1)], &[(1, 0)]]);
        assert_eq!(dual_set(&g0), regions(&[&[(-1, 0), (0, -1)]]));
    }

    #[test]
    fn type_o_is_self_dual() {
        let g0 = regions(&[&[(0, 0)]]);
        assert_eq!(dual_set(&g0), g0);
    }

    // The G0-types of Table 2 in [SS], and T2,3#T3,5.
    #[test]
    fn dual_is_an_involution() {
        let types: [&[&[(isize, isize)]]; 7] = [
            &[&[(0, 1)], &[(1, 0)]],
            &[&[(0, 2)], &[(1, 1)], &[(2, 0)]],
            &[&[(0, 2)], &[(0, 1), (1, 0)], &[(2, 0)]],
            &[&[(0, 3)], &[(1, 1)], &[(3, 0)]],
            &[&[(0, 3)], &[(1, 2)], &[(2, 1)], &[(3, 0)]],
            &[&[(0, 4)], &[(1, 2)], &[(2, 1)], &[(4, 0)]],
            &[&[(0, 5)], &[(1, 3)], &[(2, 2)], &[(3, 1)], &[(5, 0)]],
        ];
        for ty in types {
            let g0 = regions(ty);
            assert_eq!(dual_set(&dual_set(&g0)), g0);
        }
    }
}
