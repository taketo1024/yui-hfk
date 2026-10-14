use std::cmp::Ordering;

#[cfg(feature = "multithread")]
use rayon::prelude::*;

use super::f2mat::F2Mat;
use super::rng::SplitMix64;

/// `q·b` for a random unit upper-triangular `q`: each row `i` gets one random later row
/// `partner[i]` added. `q` is invertible, so `ker(q·b) = ker b`, and `(q·b)ᵀ(q·b) = bᵀ qᵀq b`
/// pairs the rows of `b` by a random form, under which repeated rows no longer cancel.
/// Applied to `bᵀ`, it mixes the columns of `b` instead (see `unmix`).
pub(super) fn mix_rows(b: &F2Mat, seed: u64) -> (F2Mat, Vec<Option<usize>>) {
    let mut rng = SplitMix64(seed ^ MIX_SEED);
    let m = b.n_rows;
    let partner: Vec<Option<usize>> = (0..m).map(|i| {
        (i + 1 < m).then(|| i + 1 + (rng.next() % (m - i - 1) as u64) as usize)
    }).collect();

    let mix = |i: usize| match partner[i] {
        Some(j) => xor_sorted(b.row(i), b.row(j)),
        None => b.row(i).to_vec(),
    };

    cfg_if::cfg_if! {
        if #[cfg(feature = "multithread")] {
            let rows: Vec<Vec<u32>> = (0..m).into_par_iter().map(mix).collect();
        } else {
            let rows: Vec<Vec<u32>> = (0..m).map(mix).collect();
        }
    }

    (F2Mat::from_rows(b.n_cols, rows), partner)
}

const MIX_SEED: u64 = 0x6d69_785f_726f_7773;

/// If column `j` of `b·r` is column `j` plus column `partner[j]` of `b`, then
/// `b·r·w = 0` iff `b·(r·w) = 0`, where `(r·w)[k] = w[k] + Σ { w[j] : partner[j] = k }`.
pub(super) fn unmix(w: &[u64], partner: &[Option<usize>]) -> Vec<u64> {
    let mut u = w.to_vec();
    for (j, &k) in partner.iter().enumerate() {
        if let Some(k) = k {
            u[k] ^= w[j];
        }
    }
    u
}

/// The symmetric difference of two increasing index lists, i.e. the sum of two 0/1 rows.
fn xor_sorted(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_utils::{SHAPES, dense, rand_block, rand_f2mat};

    #[test]
    fn xor_sorted_is_symmetric_difference() {
        let mut rng = SplitMix64(10);
        for _ in 0..100 {
            let a: Vec<u32> = (0..50).filter(|_| rng.next() & 1 == 1).collect();
            let b: Vec<u32> = (0..50).filter(|_| rng.next() & 1 == 1).collect();
            let expected: Vec<u32> = (0..50).filter(|x| a.contains(x) != b.contains(x)).collect();
            assert_eq!(xor_sorted(&a, &b), expected);
        }
        assert_eq!(xor_sorted(&[1, 2, 3], &[1, 2, 3]), Vec::<u32>::new());
        assert_eq!(xor_sorted(&[], &[4, 5]), vec![4, 5]);
    }

    #[test]
    fn mix_rows_adds_partner_rows() {
        let mut rng = SplitMix64(11);
        for (m, n) in SHAPES {
            let b = rand_f2mat(&mut rng, (m, n), 0.05);
            let (qb, partner) = mix_rows(&b, rng.next());
            let (d, dq) = (dense(&b), dense(&qb));

            assert_eq!(partner[m - 1], None);
            for i in 0..m {
                let expected: Vec<bool> = match partner[i] {
                    Some(j) => {
                        assert!(i < j && j < m);
                        d[i].iter().zip(&d[j]).map(|(x, y)| x ^ y).collect()
                    }
                    None => d[i].clone(),
                };
                assert_eq!(dq[i], expected);
            }
        }
    }

    // Mixing the rows of bᵀ gives b·r, and (b·r)·w = b·unmix(w).
    #[test]
    fn unmix_undoes_column_mixing() {
        let mut rng = SplitMix64(12);
        for (m, n) in SHAPES {
            let b = rand_f2mat(&mut rng, (m, n), 0.05);
            let (rt, partner) = mix_rows(&b.transpose(), rng.next());
            let w = rand_block(&mut rng, n);
            assert_eq!(rt.transpose().mul(&w), b.mul(&unmix(&w, &partner)));
        }
    }
}
