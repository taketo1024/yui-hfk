use log::debug;
use num_traits::{One, Zero};
use yui_core::num::FF2;

use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};

use super::f2mat::F2Mat;
use super::kernel::kernel_block;

/// Finds `w` with `aᵀw = 0` and `yᵀw = 1`, which proves `y ∉ im a`.
/// Tries up to `rounds` random starts.
pub fn find_ng_certificate(a: &SpMat<FF2>, y: &SpVec<FF2>, rounds: usize) -> Option<SpVec<FF2>> {
    let (m, n) = a.shape();
    assert_eq!(y.dim(), m);

    if y.is_zero() {
        return None;
    }

    // kernel vectors of b = aᵀ (an n × m matrix)
    let b = F2Mat::from_entries(m, n, entries_of(a)).transpose();
    let y_idx = support(y);

    (0..rounds).find_map(|round| {
        let w = kernel_block(&b, SEED + round as u64)?;
        let hits = y_idx.iter().fold(0, |acc, &i| acc ^ w[i]); // bit c = yᵀw_c
        debug!("lanczos: round {round}, certificates: {}", hits.count_ones());
        (hits != 0).then(|| column(&w, hits.trailing_zeros()))
    })
}

/// Finds `x` with `a·x = y`. Tries up to `rounds` random starts.
pub fn find_solution(a: &SpMat<FF2>, y: &SpVec<FF2>, rounds: usize) -> Option<SpVec<FF2>> {
    let (m, n) = a.shape();
    assert_eq!(y.dim(), m);

    if y.is_zero() {
        return Some(SpVec::zero(n));
    }

    // kernel vectors (x, t) of b = [a | y]; those with t = 1 give a·x = y
    let y_idx = support(y);
    let b = F2Mat::from_entries(m, n + 1, entries_of(a).chain(y_idx.iter().map(|&i| (i, n))));

    (0..rounds).find_map(|round| {
        let w = kernel_block(&b, SEED + round as u64)?;
        let hits = w[n]; // bit c = last entry of w_c
        debug!("lanczos: round {round}, solutions: {}", hits.count_ones());
        (hits != 0).then(|| column(&w[..n], hits.trailing_zeros()))
    })
}

const SEED: u64 = 0x5eed_1a9c_2b3d_4e5f;

fn entries_of(a: &SpMat<FF2>) -> impl Iterator<Item = (usize, usize)> + '_ {
    a.iter_nz().filter(|(_, _, v)| !v.is_zero()).map(|(i, j, _)| (i, j))
}

fn support(y: &SpVec<FF2>) -> Vec<usize> {
    y.iter_nz().filter(|(_, v)| !v.is_zero()).map(|(i, _)| i).collect()
}

/// Column `c` of the block `w`, as a vector.
fn column(w: &[u64], c: u32) -> SpVec<FF2> {
    let entries = w.iter().enumerate().filter(|(_, r)| (*r >> c) & 1 == 1).map(|(i, _)| (i, FF2::one()));
    SpVec::from_sorted_entries(w.len(), entries)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::rng::SplitMix64;
    use super::super::test_utils::SHAPES;
    use yui_matrix::sparse::pluq::solve_pluq;

    fn rand_mat(rng: &mut SplitMix64, (m, n): (usize, usize), density: f64) -> SpMat<FF2> {
        let threshold = (density * u64::MAX as f64) as u64;
        let entries = (0..m).flat_map(|i| (0..n).map(move |j| (i, j)))
            .filter(|_| rng.next() < threshold)
            .map(|(i, j)| (i, j, FF2::one()));
        SpMat::from_entries((m, n), entries)
    }

    fn rand_vec(rng: &mut SplitMix64, n: usize) -> SpVec<FF2> {
        SpVec::from_entries(n, (0..n).filter(|_| rng.next() & 1 == 1).map(|i| (i, FF2::one())))
    }

    fn dot(x: &SpVec<FF2>, y: &SpVec<FF2>) -> FF2 {
        x.clone().into_dense().into_iter().zip(y.clone().into_dense()).fold(FF2::zero(), |acc, (a, b)| acc + a * b)
    }

    #[test]
    fn ng_certificate_matches_pluq() {
        let mut rng = SplitMix64(1);
        let mut n_ng = 0;
        for (m, n) in SHAPES {
            for _ in 0..5 {
                let a = rand_mat(&mut rng, (m, n), 0.01);
                let y = rand_vec(&mut rng, m);
                let cert = find_ng_certificate(&a, &y, 5);

                if solve_pluq(&a, &y).is_some() {
                    assert!(cert.is_none());
                } else {
                    let w = cert.expect("no certificate found");
                    assert!((&a.transpose() * &w).is_zero());
                    assert_eq!(dot(&y, &w), FF2::one());
                    n_ng += 1;
                }
            }
        }
        assert!(n_ng > 0, "no unsolvable case was generated");
    }

    #[test]
    fn no_certificate_in_image() {
        let mut rng = SplitMix64(2);
        for (m, n) in SHAPES {
            for _ in 0..5 {
                let a = rand_mat(&mut rng, (m, n), 0.01);
                let y = &a * &rand_vec(&mut rng, n);
                assert!(find_ng_certificate(&a, &y, 5).is_none());
            }
        }
    }

    #[test]
    fn solution_matches_pluq() {
        let mut rng = SplitMix64(3);
        for (m, n) in SHAPES {
            for _ in 0..5 {
                let a = rand_mat(&mut rng, (m, n), 0.01);
                let y = &a * &rand_vec(&mut rng, n);
                let x = find_solution(&a, &y, 5).expect("no solution found");
                assert_eq!(&a * &x, y);
            }
        }
    }

    #[test]
    fn no_solution_outside_image() {
        let mut rng = SplitMix64(4);
        for (m, n) in SHAPES {
            for _ in 0..5 {
                let a = rand_mat(&mut rng, (m, n), 0.01);
                let y = rand_vec(&mut rng, m);
                if solve_pluq(&a, &y).is_none() {
                    assert!(find_solution(&a, &y, 5).is_none());
                }
            }
        }
    }

    // a = [a0 | a0]: every column appears twice, so a·aᵀ = 0 over 𝔽₂.
    #[test]
    fn duplicate_columns() {
        let mut rng = SplitMix64(5);
        for (m, n) in SHAPES {
            for _ in 0..5 {
                let a0 = rand_mat(&mut rng, (m, n), 0.01);
                let entries = a0.iter_nz().flat_map(|(i, j, _)| [(i, j, FF2::one()), (i, j + n, FF2::one())]);
                let a = SpMat::from_entries((m, 2 * n), entries);

                let y = &a * &rand_vec(&mut rng, 2 * n);
                let x = find_solution(&a, &y, 5).expect("no solution found");
                assert_eq!(&a * &x, y);

                let y = rand_vec(&mut rng, m);
                if solve_pluq(&a, &y).is_none() {
                    let w = find_ng_certificate(&a, &y, 5).expect("no certificate found");
                    assert!((&a.transpose() * &w).is_zero());
                    assert_eq!(dot(&y, &w), FF2::one());
                }
            }
        }
    }
}
