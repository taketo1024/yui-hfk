use std::fmt;

use itertools::Itertools;
use num_traits::{One, Zero};
use yui_core::num::Ratio;

use crate::g_zero::Region;

type Q = Ratio<i64>;

/// The Υ-invariant `Υ_K` on `[0, 2]`, a piecewise linear function given by its vertices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Upsilon {
    vertices: Vec<(Q, Q)>, // (t, Υ_K(t)), sorted by t, without collinear ones
}

impl Upsilon {
    /// `Υ_K(t) = −2 min{ s | ∃R ∈ G0(K), R ⊂ R_t(s) }` ([SS] Thm 2.20), where `R ⊂ R_t(s)`
    /// iff `(1 − t/2)·i + (t/2)·j ≤ s` for every corner `(i, j)` of `R`.
    pub fn from_g0(g0: &[Region]) -> Self {
        assert!(!g0.is_empty(), "G0 is non-empty");

        let f = |t: Q| -> Q {
            let min = g0.iter()
                .map(|r| r.corners().iter().map(|c| line(c.i, c.j, t)).max().unwrap())
                .min()
                .unwrap();
            -Q::from(2) * min
        };

        // note: Υ_K is linear between the intersections of the lines of the corners.
        let lines = g0.iter().flat_map(|r| r.corners()).map(|c| (c.i, c.j)).unique().collect_vec();
        let ts = lines.iter().array_combinations()
            .filter_map(|[&(i0, j0), &(i1, j1)]| intersection(i0, j0, i1, j1))
            .chain([Q::zero(), Q::from(2)])
            .filter(|t| &Q::zero() <= t && t <= &Q::from(2))
            .sorted()
            .dedup()
            .collect_vec();

        let vertices = ts.into_iter().map(|t| (t, f(t))).collect_vec();
        Upsilon { vertices: drop_collinear(vertices) }
    }

    pub fn vertices(&self) -> &[(Q, Q)] {
        &self.vertices
    }

    /// `Υ_K(t)` for `t ∈ [0, 2]`.
    pub fn value(&self, t: Q) -> Q {
        assert!(Q::zero() <= t && t <= Q::from(2), "t = {t} is not in [0, 2]");
        let ((t0, y0), (t1, y1)) = self.vertices.iter()
            .tuple_windows()
            .find(|(_, (t1, _))| t <= *t1)
            .unwrap();
        *y0 + (*y1 - y0) * (t - t0) / (*t1 - t0)
    }
}

impl fmt::Display for Upsilon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let vs = self.vertices.iter().map(|(t, y)| format!("({t}, {y})")).join(", ");
        write!(f, "PL[{vs}]")
    }
}

// The line of the corner (i, j): t ↦ (1 − t/2)·i + (t/2)·j.
fn line(i: isize, j: isize, t: Q) -> Q {
    let half = t / Q::from(2);
    (Q::one() - half) * Q::from(i as i64) + half * Q::from(j as i64)
}

// The t where the lines of (i0, j0) and (i1, j1) meet, unless they are parallel.
fn intersection(i0: isize, j0: isize, i1: isize, j1: isize) -> Option<Q> {
    let slope = ((j0 - i0) - (j1 - i1)) as i64;
    (slope != 0).then(|| Q::new(2 * (i1 - i0) as i64, slope))
}

fn drop_collinear(vs: Vec<(Q, Q)>) -> Vec<(Q, Q)> {
    let collinear = |(t0, y0): &(Q, Q), (t1, y1): &(Q, Q), (t2, y2): &(Q, Q)| {
        (*y1 - y0) * (*t2 - t1) == (*y2 - y1) * (*t1 - t0)
    };
    let inner = vs.iter().tuple_windows()
        .filter(|(a, b, c)| !collinear(a, b, c))
        .map(|(_, b, _)| *b);
    vs.first().copied().into_iter()
        .chain(inner)
        .chain(vs.last().copied().filter(|_| vs.len() > 1))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::g_zero::{Corner, dual};

    fn g0(ty: &[&[(isize, isize)]]) -> Vec<Region> {
        ty.iter().map(|r| Region::new(r.iter().map(|&(i, j)| Corner::new(i, j)))).collect()
    }

    fn q(n: i64, d: i64) -> Q {
        Q::new(n, d)
    }

    fn pl(vs: &[((i64, i64), (i64, i64))]) -> Vec<(Q, Q)> {
        vs.iter().map(|&((tn, td), (yn, yd))| (q(tn, td), q(yn, yd))).collect()
    }

    // Table 3 in [SS] (v2), where Υ_O = 0 is PL[(0, 0), (2, 0)].
    #[test]
    fn table_3() {
        let table: [(&[&[(isize, isize)]], Vec<(Q, Q)>); 6] = [
            (&[&[(0, 0)]], pl(&[((0, 1), (0, 1)), ((2, 1), (0, 1))])),
            (&[&[(0, 1)], &[(1, 0)]], pl(&[((0, 1), (0, 1)), ((1, 1), (-1, 1)), ((2, 1), (0, 1))])),
            (&[&[(0, 2)], &[(1, 1)], &[(2, 0)]], pl(&[((0, 1), (0, 1)), ((1, 1), (-2, 1)), ((2, 1), (0, 1))])),
            (&[&[(0, 2)], &[(0, 1), (1, 0)], &[(2, 0)]], pl(&[((0, 1), (0, 1)), ((2, 3), (-4, 3)), ((1, 1), (-1, 1)), ((4, 3), (-4, 3)), ((2, 1), (0, 1))])),
            (&[&[(0, 3)], &[(1, 1)], &[(3, 0)]], pl(&[((0, 1), (0, 1)), ((2, 3), (-2, 1)), ((4, 3), (-2, 1)), ((2, 1), (0, 1))])),
            (&[&[(0, 4)], &[(1, 2)], &[(2, 1)], &[(4, 0)]], pl(&[((0, 1), (0, 1)), ((2, 3), (-8, 3)), ((1, 1), (-3, 1)), ((4, 3), (-8, 3)), ((2, 1), (0, 1))])),
        ];
        for (ty, expected) in table {
            assert_eq!(Upsilon::from_g0(&g0(ty)).vertices(), expected, "{ty:?}");
        }
    }

    // Υ_{T3,4} = −3t on [0, 2/3].
    #[test]
    fn value() {
        let u = Upsilon::from_g0(&g0(&[&[(0, 3)], &[(1, 1)], &[(3, 0)]]));
        assert_eq!(u.value(q(0, 1)), q(0, 1));
        assert_eq!(u.value(q(1, 3)), q(-1, 1));
        assert_eq!(u.value(q(1, 1)), q(-2, 1));
        assert_eq!(u.value(q(2, 1)), q(0, 1));
    }

    // Υ_{K*} = −Υ_K, with G0(K*) = dual(G0(K)).
    #[test]
    fn mirror() {
        let types: [&[&[(isize, isize)]]; 4] = [
            &[&[(0, 1)], &[(1, 0)]],
            &[&[(0, 2)], &[(0, 1), (1, 0)], &[(2, 0)]],
            &[&[(0, 3)], &[(1, 1)], &[(3, 0)]],
            &[&[(0, 5)], &[(1, 3)], &[(2, 2)], &[(3, 1)], &[(5, 0)]],
        ];
        for ty in types {
            let u = Upsilon::from_g0(&g0(ty));
            let u_mirror = Upsilon::from_g0(&dual(&g0(ty)));
            let negated = u.vertices().iter().map(|&(t, y)| (t, -y)).collect_vec();
            assert_eq!(u_mirror.vertices(), negated);
        }
    }
}
