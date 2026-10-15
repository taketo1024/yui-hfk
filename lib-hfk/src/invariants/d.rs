use log::info;
use num_integer::Integer;
use num_traits::Zero;
use yui_core::num::Ratio;

use crate::g_zero::{Candidates, Region};
use crate::g_zero::calc::G0CalcError;
use crate::{G0Calc, G0CalcConfig, GridDiagram, genus};

type Q = Ratio<i64>;

pub fn d(gd: &GridDiagram, p: i64, q: i64) -> Result<Vec<Q>, G0CalcError> {
    d_with_config(gd, p, q, None, None, G0CalcConfig::default())
}

/// `[d(S³_{p/q}(K), i) | 0 ≤ i < |p|]` for coprime `p ≠ 0` and `q > 0`. For `p < 0`, it is
/// `−d(S³_{−p/q}(K*), i)` by `S³_{p/q}(K) ≅ −S³_{−p/q}(K*)`, with the same labels `i`.
pub fn d_with_config(gd: &GridDiagram, p: i64, q: i64, g: Option<usize>, g4: Option<usize>, config: G0CalcConfig) -> Result<Vec<Q>, G0CalcError> {
    assert!(p != 0, "0-surgery is not supported");
    if p > 0 {
        let vs = v_seq_with_config(gd, g, g4, config)?;
        Ok(d_from_v(&vs, p, q))
    } else {
        let vs = v_seq_with_config(&gd.rotate90(), g, g4, config)?;
        Ok(d_from_v(&vs, -p, q).into_iter().map(|d| -d).collect())
    }
}

/// `d(S³_{p/q}(K), i) = d(S³_{p/q}(O), i) − 2·max{ V_⌊i/q⌋(K), V_⌊(p+q−1−i)/q⌋(K) }` for
/// `0 ≤ i < p` ([SS] Prop 2.19, [NW15, Prop 1.6]), where `vs = [V_0(K), V_1(K), …]` and
/// `V_k(K) = 0` beyond.
pub fn d_from_v(vs: &[usize], p: i64, q: i64) -> Vec<Q> {
    assert!(p > 0 && q > 0 && p.gcd(&q) == 1, "p/q = {p}/{q} must be positive and reduced");
    let v = |k: i64| vs.get(k as usize).copied().unwrap_or(0) as i64;
    (0..p)
        .map(|i| d_lens(p, q, i) - Q::from(2 * v(i / q).max(v((p + q - 1 - i) / q))))
        .collect()
}

/// `[V_0(K), …, V_{g4−1}(K)]`, as `V_k(K) = 0` for `k ≥ g4`. With `V_k(K) = min{ m ≥ 0 |
/// R(m, k+m) is realizable }` ([SS] Thm 2.13, 2.20) and `V_k(K) ≤ ⌈(g4 − k)/2⌉` (Rasmussen),
/// these simple regions are tested all at once.
pub fn v_seq_with_config(gd: &GridDiagram, g: Option<usize>, g4: Option<usize>, config: G0CalcConfig) -> Result<Vec<usize>, G0CalcError> {
    info!("compute V_k ...");

    let g4 = g4.unwrap_or_else(|| g.unwrap_or_else(|| genus(gd)));
    let bound = |k: usize| (g4 - k).div_ceil(2);
    let region = |k: usize, m: usize| Region::simple(m as isize, (k + m) as isize);

    let mut candidates = Candidates::new();
    (0..g4).flat_map(|k| (0..=bound(k)).map(move |m| (k, m)))
        .for_each(|(k, m)| candidates.add(region(k, m)));

    let mut calc = G0Calc::init_with_candidates(gd, candidates, config);
    calc.run();

    let realizable = |r: &Region| calc.accepted().iter().any(|a| r.includes(a));
    let unrealizable = |r: &Region| calc.rejected().iter().any(|s| s.includes(r));

    (0..g4).map(|k| {
        let v = (0..=bound(k)).find(|&m| realizable(&region(k, m)))
            .unwrap_or_else(|| panic!("Invalid state: no R(m, {k}+m) with m ≤ {} is realizable", bound(k)));
        if (0..v).all(|m| unrealizable(&region(k, m))) {
            Ok(v)
        } else {
            Err(G0CalcError::Resigned)
        }
    }).collect()
}

/// `V_k(K) = min{ m ≥ 0 | ∃R ∈ G0(K), R ⊂ R(m, k+m) }` ([SS] Thm 2.20), where `R ⊂ R(m, k+m)`
/// iff `i ≤ m` and `j ≤ k + m` for every corner `(i, j)` of `R`.
pub fn v_from_g0(g0: &[Region], k: usize) -> usize {
    let k = k as isize;
    let m = g0.iter()
        .map(|r| r.corners().iter().map(|c| c.i.max(c.j - k)).max().unwrap_or(isize::MIN))
        .min()
        .expect("G0 is non-empty");
    m.max(0) as usize
}

/// `d(S³_{p/q}(O), i) = −1/4 + (2i + 1 − p − q)²/(4pq) − d(S³_{q/r}(O), j)`, where `r` and `j`
/// are the residues of `p` and `i` mod `q` ([OS03a, Prop 4.8], in the convention of [NW15]).
fn d_lens(p: i64, q: i64, i: i64) -> Q {
    if q == 0 {
        return Q::zero(); // S³_{1/0}(O) = S³
    }
    let a = 2 * i + 1 - p - q;
    Q::new(-1, 4) + Q::new(a * a, 4 * p * q) - d_lens(q, p % q, i % q)
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;
    use crate::g_zero::{Corner, dual};

    fn g0(ty: &[&[(isize, isize)]]) -> Vec<Region> {
        ty.iter().map(|r| Region::new(r.iter().map(|&(i, j)| Corner::new(i, j)))).collect()
    }

    fn q(n: i64, d: i64) -> Q {
        Q::new(n, d)
    }

    // d(S³_p(O), i) = ((2i − p)² − p)/(4p).
    #[test]
    fn lens_integral() {
        for p in 1..=6 {
            let expected = (0..p).map(|i| q((2 * i - p) * (2 * i - p) - p, 4 * p)).collect_vec();
            assert_eq!(d_from_v(&[], p, 1), expected, "p = {p}");
        }
    }

    // T2,3: V_0 = 1, so d(S³_1) = −2 and d(S³_2) = (1/4 − 2, −1/4).
    #[test]
    fn t2_3() {
        assert_eq!(d_from_v(&[1], 1, 1), [q(-2, 1)]);
        assert_eq!(d_from_v(&[1], 2, 1), [q(-7, 4), q(-1, 4)]);
    }

    // V_k(K) and V_k(K*) of Table 2 in [SS] (v2), with V_k = 0 beyond.
    #[test]
    fn table_2() {
        let table: [(&[&[(isize, isize)]], &[usize]); 6] = [
            (&[&[(0, 0)]], &[0]),
            (&[&[(0, 1)], &[(1, 0)]], &[1, 0]),
            (&[&[(0, 2)], &[(1, 1)], &[(2, 0)]], &[1, 1, 0]),
            (&[&[(0, 2)], &[(0, 1), (1, 0)], &[(2, 0)]], &[1, 1, 0]),
            (&[&[(0, 3)], &[(1, 1)], &[(3, 0)]], &[1, 1, 1, 0]),
            (&[&[(0, 4)], &[(1, 2)], &[(2, 1)], &[(4, 0)]], &[2, 1, 1, 1, 0]),
        ];
        for (ty, vs) in table {
            let n = vs.len() + 1;
            let expected = vs.iter().copied().chain([0]).collect_vec();
            assert_eq!((0..n).map(|k| v_from_g0(&g0(ty), k)).collect_vec(), expected, "{ty:?}");
            assert!((0..n).all(|k| v_from_g0(&dual(&g0(ty)), k) == 0), "{ty:?}: mirror");
        }
    }

    // note: Table 2 leaves V_k(T2,7) blank. They agree with the torsion coefficients.
    #[test]
    fn t2_7() {
        let t2_7 = g0(&[&[(0, 3)], &[(1, 2)], &[(2, 1)], &[(3, 0)]]);
        assert_eq!((0..4).map(|k| v_from_g0(&t2_7, k)).collect_vec(), [2, 1, 1, 0]);
    }

    // V_k(K*) = min{ m ≥ 0 | ∀R ∈ G0(K), (−m, −k − m) ∈ R } ([SS] Prop 2.21).
    #[test]
    fn v_of_mirror() {
        let t3_4 = g0(&[&[(0, 3)], &[(1, 1)], &[(3, 0)]]);
        let mirror = dual(&t3_4);
        for k in 0..4 {
            let m = (0..).find(|&m: &isize| t3_4.iter().all(|r| r.contains((-m, -(k as isize) - m)))).unwrap();
            assert_eq!(v_from_g0(&mirror, k), m as usize, "k = {k}");
        }
    }

    // T(2,5) = 5_1, with G0-type T2,5: V = (1, 1), and V = (0, 0) for the mirror.
    #[test]
    fn v_seq_5_1() {
        let gd = GridDiagram::test_data("5_1").unwrap();
        let config = G0CalcConfig::default();
        assert_eq!(v_seq_with_config(&gd, Some(2), Some(2), config), Ok(vec![1, 1]));
        assert_eq!(v_seq_with_config(&gd.rotate90(), Some(2), Some(2), config), Ok(vec![0, 0]));

        let t2_5 = g0(&[&[(0, 2)], &[(1, 1)], &[(2, 0)]]);
        assert_eq!((0..2).map(|k| v_from_g0(&t2_5, k)).collect_vec(), [1, 1]);
    }

    #[test]
    fn d_trefoil() {
        assert_eq!(d(&GridDiagram::trefoil(), 1, 1), Ok(vec![q(-2, 1)]));
    }

    #[test]
    fn d_5_1() {
        let gd = GridDiagram::test_data("5_1").unwrap();
        let config = G0CalcConfig::default();
        assert_eq!(d_with_config(&gd, 1, 1, Some(2), Some(2), config), Ok(vec![q(-2, 1)]));
        assert_eq!(d_with_config(&gd, -1, 1, Some(2), Some(2), config), Ok(vec![q(0, 1)]));
    }
}
