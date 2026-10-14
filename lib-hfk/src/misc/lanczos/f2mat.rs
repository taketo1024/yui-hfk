use std::iter::once;

#[cfg(feature = "multithread")]
use rayon::prelude::*;

/// A sparse 0/1 matrix in CSR form.
pub(super) struct F2Mat {
    pub(super) n_rows: usize,
    pub(super) n_cols: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<u32>,
}

impl F2Mat {
    pub(super) fn from_entries(n_rows: usize, n_cols: usize, entries: impl IntoIterator<Item = (usize, usize)>) -> Self {
        assert!(n_cols <= u32::MAX as usize);
        let entries: Vec<(usize, usize)> = entries.into_iter().collect();

        let mut row_ptr = vec![0; n_rows + 1];
        for &(i, _) in &entries {
            row_ptr[i + 1] += 1;
        }
        for i in 0..n_rows {
            row_ptr[i + 1] += row_ptr[i];
        }

        let mut next = row_ptr.clone();
        let mut col_idx = vec![0; entries.len()];
        for &(i, j) in &entries {
            col_idx[next[i]] = j as u32;
            next[i] += 1;
        }

        Self { n_rows, n_cols, row_ptr, col_idx }
    }

    pub(super) fn from_rows(n_cols: usize, rows: Vec<Vec<u32>>) -> Self {
        let row_ptr = once(0).chain(rows.iter().scan(0, |acc, r| {
            *acc += r.len();
            Some(*acc)
        })).collect();
        Self { n_rows: rows.len(), n_cols, row_ptr, col_idx: rows.concat() }
    }

    pub(super) fn transpose(&self) -> Self {
        let entries = (0..self.n_rows).flat_map(|i| self.row(i).iter().map(move |&j| (j as usize, i)));
        Self::from_entries(self.n_cols, self.n_rows, entries)
    }

    pub(super) fn row(&self, i: usize) -> &[u32] {
        &self.col_idx[self.row_ptr[i]..self.row_ptr[i + 1]]
    }

    /// The block product `self·v`, for a block `v` with one `u64` per column.
    pub(super) fn mul(&self, v: &[u64]) -> Vec<u64> {
        assert_eq!(v.len(), self.n_cols);
        let f = |i: usize| self.row(i).iter().fold(0, |acc, &j| acc ^ v[j as usize]);

        cfg_if::cfg_if! {
            if #[cfg(feature = "multithread")] {
                (0..self.n_rows).into_par_iter().map(f).collect()
            } else {
                (0..self.n_rows).map(f).collect()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::mat64::inner;
    use super::super::rng::SplitMix64;
    use super::super::test_utils::{SHAPES, dense, rand_block, rand_f2mat};

    #[test]
    fn f2mat_from_entries_and_rows_agree() {
        let mut rng = SplitMix64(6);
        let a = rand_f2mat(&mut rng, (300, 200), 0.05);
        let b = F2Mat::from_rows(a.n_cols, (0..a.n_rows).map(|i| a.row(i).to_vec()).collect());
        assert_eq!((b.n_rows, b.n_cols, &b.row_ptr, &b.col_idx), (a.n_rows, a.n_cols, &a.row_ptr, &a.col_idx));

        // empty rows
        let c = F2Mat::from_entries(3, 4, [(0, 1), (2, 0), (2, 3)]);
        assert_eq!(dense(&c), [
            vec![false, true, false, false],
            vec![false, false, false, false],
            vec![true, false, false, true],
        ]);
    }

    #[test]
    fn f2mat_transpose() {
        let mut rng = SplitMix64(7);
        for (m, n) in SHAPES {
            let a = rand_f2mat(&mut rng, (m, n), 0.05);
            let at = a.transpose();
            let (d, dt) = (dense(&a), dense(&at));

            assert_eq!((at.n_rows, at.n_cols), (n, m));
            assert!((0..m).all(|i| (0..n).all(|j| d[i][j] == dt[j][i])));
            assert!((0..n).all(|j| at.row(j).windows(2).all(|w| w[0] < w[1])), "rows are sorted");
            assert_eq!(dense(&at.transpose()), d);
        }
    }

    #[test]
    fn f2mat_mul() {
        let mut rng = SplitMix64(8);
        for (m, n) in SHAPES {
            let a = rand_f2mat(&mut rng, (m, n), 0.05);
            let v = rand_block(&mut rng, n);
            let expected: Vec<u64> = dense(&a).iter().map(|row| {
                row.iter().zip(&v).filter(|(e, _)| **e).fold(0, |acc, (_, x)| acc ^ x)
            }).collect();
            assert_eq!(a.mul(&v), expected);
        }
    }

    // uᵀ(a·v) = (aᵀu)ᵀv, which ties transpose, mul and inner together.
    #[test]
    fn f2mat_mul_transpose_adjoint() {
        let mut rng = SplitMix64(9);
        for (m, n) in SHAPES {
            let a = rand_f2mat(&mut rng, (m, n), 0.05);
            let (u, v) = (rand_block(&mut rng, m), rand_block(&mut rng, n));
            assert_eq!(inner(&u, &a.mul(&v)), inner(&a.transpose().mul(&u), &v));
        }
    }
}
