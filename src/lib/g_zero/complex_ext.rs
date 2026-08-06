use std::ops::RangeInclusive;

use num_integer::binomial;
use yui_core::num::FF2;
use yui_homology::Summand;

use crate::grid_complex::{GCFiltered, GCGen, GCInfGen, GridComplex};
use super::Region;

// An inflated generator U^a x has
//
//   maslov = M(x) - 2e,   alex = A(x) - e,   alg = -a₀,      e = Σᵢ aᵢ.
//
// So alg does *not* determine e: a point (alg, alex) = (i, j) at maslov degree k
// draws from every total degree e ≥ -i, taking raw gens at (k + 2e, j + e) and
// monomials with U₀-exponent exactly -i.
impl GridComplex<GCFiltered> {
    // Returns the (alg_deg, alex_deg) coordinates outside r that can host inf_gens at maslov_deg k.
    pub(crate) fn complementary_points(&self, k: isize, r: &Region) -> Vec<(isize, isize)> {
        let n = self.grid_number() as usize;

        let mut pts: Vec<(isize, isize)> = self.gens.support()
            .filter_map(|&(m, a)| {
                let d = m - k; // = 2e
                (d >= 0 && d % 2 == 0).then_some((a, d / 2))
            })
            .flat_map(|(a, e)| {
                let j = a - e;
                (-e ..= 0).map(move |i| (i, j, e))
            })
            .filter(|&(i, j, e)| {
                n_monomials(n, (e + i) as usize) > 0 && !r.contains((i, j))
            })
            .map(|(i, j, _)| (i, j))
            .collect();

        pts.sort_unstable();
        pts.dedup();
        pts
    }

    // Count GCInfGens of maslov_deg = k at the given (alg, alex) points without allocating them.
    pub(crate) fn count_inf_gens_in(&self, k: isize, pts: &[(isize, isize)]) -> usize {
        let n = self.grid_number() as usize;
        self.contributions(k, pts).map(|(i, j, e)| {
            let n_bases = self.gens.raw_gens(k + 2 * e, j + e).len();
            let n_multi = n_monomials(n, (e + i) as usize);
            n_bases * n_multi
        }).sum()
    }

    // Collect all GCInfGens of maslov_deg = k at the given (alg, alex) points.
    pub(crate) fn inf_gens_in(&self, k: isize, pts: &[(isize, isize)]) -> Summand<GCInfGen, FF2> {
        let gens = self.contributions(k, pts).flat_map(|(i, j, e)| {
            self.inf_gens_e_at(k, j, e as usize, (-i) as usize)
        });
        Summand::from_raw_generators(gens)
    }

    // Count base GCGens at the given points without allocating them.
    pub(crate) fn count_raw_gens_in(&self, k: isize, pts: &[(isize, isize)]) -> usize {
        self.raw_bidegrees_in(k, pts).into_iter().map(|(m, a)| {
            self.gens.raw_gens(m, a).len()
        }).sum()
    }

    // Collect all base GCGens whose (maslov, alex) corresponds to one of the given
    // inflated (alg, alex) points at maslov_deg = k.
    pub(crate) fn raw_gens_in(&self, k: isize, pts: &[(isize, isize)]) -> Summand<GCGen, FF2> {
        let gens = self.raw_bidegrees_in(k, pts).into_iter().flat_map(|(m, a)| {
            self.gens.raw_gens(m, a).iter().cloned()
        });
        Summand::from_raw_generators(gens)
    }

    // ---- internal ----

    // The (point, total-degree) triples (i, j, e) that carry generators of maslov
    // degree k. Distinct triples yield disjoint generator sets, since e is itself
    // a function of the generator.
    fn contributions<'a>(&'a self, k: isize, pts: &'a [(isize, isize)])
        -> impl Iterator<Item = (isize, isize, isize)> + 'a
    {
        let n = self.grid_number() as usize;
        pts.iter().flat_map(move |&(i, j)| {
            self.e_range(k, i)
                .filter(move |&e| n_monomials(n, (e + i) as usize) > 0)
                .map(move |e| (i, j, e))
        })
    }

    // Total U-degrees that can host a generator of maslov degree k with U₀-exponent -i:
    // needs e ≥ -i, and a raw gen at maslov k + 2e.
    fn e_range(&self, k: isize, i: isize) -> RangeInclusive<isize> {
        let m_max = *self.gens.maslov_deg_range().end();
        (-i).max(0) ..= (m_max - k).div_euclid(2)
    }

    // The distinct raw (maslov, alex) bidegrees underlying the inflated generators.
    fn raw_bidegrees_in(&self, k: isize, pts: &[(isize, isize)]) -> Vec<(isize, isize)> {
        let mut v: Vec<(isize, isize)> = self.contributions(k, pts)
            .map(|(_, j, e)| (k + 2 * e, j + e))
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    }
}

// Number of monomials in U₁ … U_{n-1} of total degree m.
fn n_monomials(n: usize, m: usize) -> usize {
    if n < 2 {
        usize::from(m == 0)
    } else {
        binomial(n - 2 + m, m)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use yui_core::poly::MultiVar;

    use super::*;
    use crate::{GridComplex, GridDiagram};
    use super::super::{Corner, Region};

    // Independent enumeration of every inflated generator of maslov degree k
    // lying outside r: all raw gens against all monomials, no indexing tricks.
    // Used to cross-check `complementary_points` + `inf_gens_in`.
    fn brute_force(cx: &GridComplex<GCFiltered>, k: isize, r: &Region) -> Vec<GCInfGen> {
        let n = cx.grid_number() as usize;
        let m_max = *cx.gens.maslov_deg_range().end();

        if m_max < k {
            return vec![];
        }

        let e_max = ((m_max - k) / 2) as usize;
        let mut out: Vec<GCInfGen> = (0..=e_max).flat_map(|e| {
            cx.gens.iter()
                .filter(move |x| x.maslov_deg == k + 2 * e as isize)
                .flat_map(move |x| {
                    (0..n).combinations_with_replacement(e).map(move |combo| {
                        let u = MultiVar::from_iter(
                            (0..n).map(|v| (v, combo.iter().filter(|&&w| w == v).count()))
                        );
                        GCInfGen::new(u, x.clone())
                    })
                })
        })
        .filter(|g| !r.contains_gen(g))
        .collect();

        out.sort();
        out
    }

    fn test_regions() -> Vec<Region> {
        vec![
            Region::simple(0, 0),
            Region::simple(-1, -1),
            Region::simple(0, 1),
            Region::simple(-1, 2),
            Region::new([Corner::new(-1, isize::MAX), Corner::new(0, 0)]),
            Region::new([Corner::new(-2, 1), Corner::new(0, -1)]),
        ]
    }

    fn test_complexes() -> Vec<GridComplex<GCFiltered>> {
        [GridDiagram::unknot(), GridDiagram::trefoil(), GridDiagram::figure8()]
            .into_iter()
            .map(GridComplex::<GCFiltered>::new)
            .collect()
    }

    // ---- cross-checks against brute force ----

    #[test]
    fn test_inf_gens_in_matches_brute_force() {
        for cx in test_complexes() {
            let m_range = cx.gens.maslov_deg_range();
            for k in (*m_range.start() - 2) ..= *m_range.end() {
                for r in test_regions() {
                    let pts = cx.complementary_points(k, &r);
                    let mut actual: Vec<GCInfGen> = cx.inf_gens_in(k, &pts)
                        .raw_generators().iter().cloned().collect();
                    actual.sort();

                    assert_eq!(actual, brute_force(&cx, k, &r), "k = {k}, r = {r}");
                }
            }
        }
    }

    #[test]
    fn test_count_inf_gens_in_matches_inf_gens_in() {
        for cx in test_complexes() {
            let m_range = cx.gens.maslov_deg_range();
            for k in (*m_range.start() - 2) ..= *m_range.end() {
                for r in test_regions() {
                    let pts = cx.complementary_points(k, &r);
                    assert_eq!(
                        cx.count_inf_gens_in(k, &pts),
                        cx.inf_gens_in(k, &pts).n_generators(),
                        "k = {k}, r = {r}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_count_raw_gens_in_matches_raw_gens_in() {
        for cx in test_complexes() {
            let m_range = cx.gens.maslov_deg_range();
            for k in (*m_range.start() - 2) ..= *m_range.end() {
                for r in test_regions() {
                    let pts = cx.complementary_points(k, &r);
                    assert_eq!(
                        cx.count_raw_gens_in(k, &pts),
                        cx.raw_gens_in(k, &pts).n_generators(),
                        "k = {k}, r = {r}"
                    );
                }
            }
        }
    }

    // Every generated element must sit at maslov degree k, at one of the
    // reported points, and outside r.
    #[test]
    fn test_inf_gens_in_are_well_placed() {
        for cx in test_complexes() {
            let m_range = cx.gens.maslov_deg_range();
            for k in (*m_range.start() - 2) ..= *m_range.end() {
                for r in test_regions() {
                    let pts = cx.complementary_points(k, &r);
                    for g in cx.inf_gens_in(k, &pts).raw_generators() {
                        assert_eq!(g.maslov_deg(), k);
                        assert!(pts.contains(&(g.alg_deg(), g.alex_deg())));
                        assert!(!r.contains_gen(g));
                    }
                }
            }
        }
    }

    #[test]
    fn test_complementary_points_are_outside_region() {
        for cx in test_complexes() {
            let m_range = cx.gens.maslov_deg_range();
            for k in (*m_range.start() - 2) ..= *m_range.end() {
                for r in test_regions() {
                    let pts = cx.complementary_points(k, &r);
                    assert!(pts.iter().all(|&p| !r.contains(p)));
                    assert!(pts.iter().all(|&(i, _)| i <= 0), "alg degrees are ≤ 0");
                    assert!(pts.iter().tuple_windows().all(|(a, b)| a < b), "sorted, no dups");
                }
            }
        }
    }

    // ---- complementary_points tests ----

    #[test]
    fn test_complementary_points_unknot() {
        // unknot: raw gens at (maslov=0, alex=0) and (maslov=-1, alex=-1), so j_max=0.
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCFiltered>::new(gd);

        // k=0: only raw gen at (0,0). R(-1,-1) does not contain (0,0) → point (0,0).
        let pts = cx.complementary_points(0, &Region::simple(-1, -1));
        assert_eq!(pts, vec![(0, 0)]);

        // k=-1: only raw gen at (-1,-1). R(-1,-1) does not contain (0,-1) → point (0,-1).
        let pts = cx.complementary_points(-1, &Region::simple(-1, -1));
        assert_eq!(pts, vec![(0, -1)]);

        // R(0,0): contains all (i,j) with i ≤ 0 and j ≤ 0 → every point is inside r.
        let pts = cx.complementary_points(0, &Region::simple(0, 0));
        assert!(pts.is_empty());
        let pts = cx.complementary_points(-1, &Region::simple(0, 0));
        assert!(pts.is_empty());
    }

    // ---- inf_gens_in tests ----

    #[test]
    fn test_inf_gens_in_unknot() {
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCFiltered>::new(gd);

        // The unknot has a single raw gen at (maslov=0, alex=0).
        // Its only inf_gen at (k=0, j=0) is the trivial one: alg=0, alex=0.
        // R(0,0) contains (0,0), so the points are empty.
        let pts = cx.complementary_points(0, &Region::simple(0, 0));
        let gen0 = cx.inf_gens_in(0, &pts);
        assert_eq!(gen0.n_generators(), 0);

        // R(-1,-1) does not contain (0,0), so we get 1 gen.
        let pts = cx.complementary_points(0, &Region::simple(-1, -1));
        let gen1 = cx.inf_gens_in(0, &pts);
        assert_eq!(gen1.n_generators(), 1);
    }

}
