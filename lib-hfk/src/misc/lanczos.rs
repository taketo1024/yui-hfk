//! Block Lanczos over `𝔽₂` (P. Montgomery, *A Block Lanczos Algorithm for Finding
//! Dependencies over GF(2)*, EUROCRYPT 1995), for sparse systems too large for PLUQ.
//!
//! [`find_ng_certificate`] looks for `w` with `aᵀw = 0` and `yᵀw = 1` (so `y ∉ im a`),
//! and [`find_solution`] for `x` with `a·x = y`. Both are randomized, but every vector
//! they return is verified, so a `Some` is always correct and a `None` only means
//! that nothing was found.

use std::array::from_fn;
use std::cmp::Ordering;
use std::iter::once;

use log::debug;
use num_traits::{One, Zero};
use yui_core::num::FF2;

use yui_matrix::MatTrait;
use yui_matrix::sparse::{SpMat, SpVec};

#[cfg(feature = "multithread")]
use rayon::prelude::*;

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

/// One run of Block Lanczos on `s = bᵀb`. Returns 64 vectors of `ker b` (some may be
/// zero) as a block with one `u64` per coordinate, or `None` on a breakdown.
fn kernel_block(b0: &F2Mat, seed: u64) -> Option<Vec<u64>> {
    // note: Lanczos runs on b = q·b0·r for random unit-triangular q, r (see `mix_rows`), so that
    // repeated rows or columns of b0 do not cancel in bᵀb. Then ker b0 = r·ker b.
    let qb0 = mix_rows(b0, seed).0;
    let (bt, partner) = mix_rows(&qb0.transpose(), !seed);
    let b = &bt.transpose();
    let bt = &bt;

    let n = b.n_cols;
    let s = |v: &[u64]| bt.mul(&b.mul(v));

    let mut rng = SplitMix64(seed);
    let y: Vec<u64> = (0..n).map(|_| rng.next()).collect();
    let v0 = s(&y); // solving s·x = s·y gives x − y ∈ ker s

    let mut v = [v0.clone(), vec![0; n], vec![0; n]]; // V_i, V_{i-1}, V_{i-2}
    let mut winv = [[0; 64]; 2];                       // W_{i-1}^inv, W_{i-2}^inv
    let mut vt_sv_prev = [0; 64];                      // V_{i-1}ᵀ s V_{i-1}
    let mut vt_s2v_prev = [0; 64];                     // V_{i-1}ᵀ s² V_{i-1}
    let mut mask_prev = u64::MAX;                      // S_{i-1}, as a column mask
    let mut x = vec![0; n];

    let max_iter = n / 60 + 100;
    let mut iter = 0;
    let mut dim = 0; // Σ |S_i|, the dimension of the explored Krylov space

    loop {
        if iter > max_iter {
            debug!("lanczos: no convergence after {iter} iterations");
            return None;
        }

        let sv = s(&v[0]);
        let vt_sv = inner(&v[0], &sv);
        if vt_sv == [0; 64] {
            break;
        }
        let vt_s2v = inner(&sv, &sv);

        // note: near the end, T_i may be too small to cover the columns skipped at the
        // previous step. We stop there; `combine` checks every vector it returns.
        let Some((mask, winv0)) = select_columns(&vt_sv, mask_prev) else {
            debug!("lanczos: column selection failed at iter {iter}, stopping");
            break;
        };
        dim += mask.count_ones() as usize;

        // the coefficients D, E, F of the recurrence (over 𝔽₂, − = +)
        let d = add_id(&mul64(&winv0, &from_fn(|i| (vt_s2v[i] & mask) ^ vt_sv[i])));
        let e = mul64(&winv[0], &vt_sv).map(|r| r & mask);
        let f = mul64(
            &mul64(&winv[1], &add_id(&mul64(&vt_sv_prev, &winv[0]))),
            &from_fn(|i| ((vt_s2v_prev[i] & mask_prev) ^ vt_sv_prev[i]) & mask),
        );

        let mut v_next: Vec<u64> = sv.iter().map(|r| r & mask).collect();
        mul_acc(&mut v_next, &v[0], &d);
        mul_acc(&mut v_next, &v[1], &e);
        mul_acc(&mut v_next, &v[2], &f);

        // X += V_i W_i^inv V_iᵀ V_0
        mul_acc(&mut x, &v[0], &mul64(&winv0, &inner(&v[0], &v0)));

        v.rotate_right(1);
        v[0] = v_next;
        winv = [winv0, winv[0]];
        vt_sv_prev = vt_sv;
        vt_s2v_prev = vt_s2v;
        mask_prev = mask;

        iter += 1;
        if iter % 1000 == 0 {
            debug!("lanczos: iter {iter} (≤ {} expected)", n / 63);
        }
    }

    debug!("lanczos: converged after {iter} iterations, explored dimension {dim}");

    let z_lo: Vec<u64> = x.iter().zip(&y).map(|(a, b)| a ^ b).collect();
    let w = unmix(&combine(b, &z_lo, &v[0])?, &partner);

    // an exact check against b0 (it holds by construction)
    let bad = b0.mul(&w).iter().fold(0, |acc, r| acc | r);
    if bad != 0 {
        debug!("lanczos: {} vectors failed the check against the unmixed matrix", bad.count_ones());
    }
    Some(w.iter().map(|r| r & !bad).collect())
}

/// `q·b` for a random unit upper-triangular `q`: each row `i` gets one random later row
/// `partner[i]` added. `q` is invertible, so `ker(q·b) = ker b`, and `(q·b)ᵀ(q·b) = bᵀ qᵀq b`
/// pairs the rows of `b` by a random form, under which repeated rows no longer cancel.
/// Applied to `bᵀ`, it mixes the columns of `b` instead (see `unmix`).
fn mix_rows(b: &F2Mat, seed: u64) -> (F2Mat, Vec<Option<usize>>) {
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
fn unmix(w: &[u64], partner: &[Option<usize>]) -> Vec<u64> {
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

/// Montgomery's choice of `S_i`: a column mask `S` and `W^inv = S (Sᵀ t S)⁻¹ Sᵀ`, where `S`
/// must contain every column not in `S_{i-1}`. `None` on failure (a breakdown).
fn select_columns(t: &Mat64, mask_prev: u64) -> Option<(u64, Mat64)> {
    let mut m: [(u64, u64); 64] = from_fn(|i| (t[i], 1 << i)); // [t | I]
    let order: Vec<usize> = (0..64).filter(|&c| (mask_prev >> c) & 1 == 0)
        .chain((0..64).filter(|&c| (mask_prev >> c) & 1 == 1))
        .collect();
    let mut mask = 0;

    for i in 0..64 {
        let c = order[i];
        let bit = 1 << c;

        if let Some(j) = (i..64).find(|&j| m[order[j]].0 & bit != 0) {
            m.swap(c, order[j]);
            eliminate(&mut m, c, |row| row.0 & bit != 0);
            mask |= bit;
        } else {
            let Some(j) = (i..64).find(|&j| m[order[j]].1 & bit != 0) else {
                debug!("lanczos: select_columns: no pivot in the inverse half");
                return None;
            };
            m.swap(c, order[j]);
            eliminate(&mut m, c, |row| row.1 & bit != 0);
            m[c] = (0, 0);
        }
    }

    if mask == 0 || mask | mask_prev != u64::MAX {
        debug!("lanczos: select_columns: |S| = {}, |S_prev| = {}, uncovered = {}",
            mask.count_ones(), mask_prev.count_ones(), (!(mask | mask_prev)).count_ones());
        return None;
    }

    Some((mask, from_fn(|i| m[i].1)))
}

fn eliminate(m: &mut [(u64, u64); 64], p: usize, hit: impl Fn(&(u64, u64)) -> bool) {
    let piv = m[p];
    for r in (0..64).filter(|&r| r != p) {
        if hit(&m[r]) {
            m[r].0 ^= piv.0;
            m[r].1 ^= piv.1;
        }
    }
}

/// Montgomery's final step: with `Z = [X − Y | V_m]`, returns `Z·u` for the `u` with `b·Z·u = 0`.
fn combine(b: &F2Mat, z_lo: &[u64], z_hi: &[u64]) -> Option<Vec<u64>> {
    let (bz_lo, bz_hi) = (b.mul(z_lo), b.mul(z_hi));
    let mut rows: Vec<u128> = bz_lo.iter().zip(&bz_hi).map(|(&l, &h)| l as u128 | (h as u128) << 64).collect();
    let tags: Vec<u128> = null_space_128(&mut rows).into_iter().take(64).collect();

    if tags.is_empty() {
        return None;
    }

    // u_lo, u_hi: the 64 × 64 halves of the matrix whose columns are the tags
    let u_half = |offset: usize| -> Mat64 {
        from_fn(|i| tags.iter().enumerate().fold(0, |acc, (t, tag)| acc | (((tag >> (i + offset)) & 1) as u64) << t))
    };
    let mut w = vec![0; z_lo.len()];
    mul_acc(&mut w, z_lo, &u_half(0));
    mul_acc(&mut w, z_hi, &u_half(64));

    // an exact check: drop any column not in ker b (none, unless something is wrong)
    let bad = b.mul(&w).iter().fold(0, |acc, r| acc | r);
    if bad != 0 {
        debug!("lanczos: {} vectors failed the kernel check", bad.count_ones());
    }

    Some(w.iter().map(|r| r & !bad).collect())
}

/// Null vectors of the tall matrix with the given rows (bit `c` = column `c`), by column elimination.
fn null_space_128(rows: &mut [u128]) -> Vec<u128> {
    let mut tags: [u128; 128] = from_fn(|c| 1 << c);
    let mut active = u128::MAX;

    for r in 0..rows.len() {
        let hit = rows[r] & active;
        if hit == 0 {
            continue;
        }

        let p = hit.trailing_zeros() as usize;
        let others = hit & !(1 << p);

        // add column p to the columns in `others`; rows above r are zero in column p.
        for row in rows[r..].iter_mut().filter(|row| (**row >> p) & 1 == 1) {
            *row ^= others;
        }
        for q in (0..128).filter(|&q| (others >> q) & 1 == 1) {
            tags[q] ^= tags[p];
        }

        active &= !(1 << p);
        if active == 0 {
            break;
        }
    }

    (0..128).filter(|&c| (active >> c) & 1 == 1).map(|c| tags[c]).collect()
}

/// A 64 × 64 matrix over `𝔽₂`: row `i` is a word, bit `j` is column `j`.
type Mat64 = [u64; 64];

fn mul64(a: &Mat64, b: &Mat64) -> Mat64 {
    from_fn(|i| (0..64).filter(|&j| (a[i] >> j) & 1 == 1).fold(0, |acc, j| acc ^ b[j]))
}

fn add_id(a: &Mat64) -> Mat64 {
    from_fn(|i| a[i] ^ (1 << i))
}

/// For blocks `a`, `b` (one `u64` per row), returns `aᵀb`.
fn inner(a: &[u64], b: &[u64]) -> Mat64 {
    assert_eq!(a.len(), b.len());
    const CHUNK: usize = 1 << 14;

    // t[k][x] = ⊕ { b[r] : byte k of a[r] is x }
    let partial = |(xs, ys): (&[u64], &[u64])| -> Tables {
        let mut t = new_tables();
        for (&x, &y) in xs.iter().zip(ys) {
            for (k, tk) in t.iter_mut().enumerate() {
                tk[((x >> (8 * k)) & 0xff) as usize] ^= y;
            }
        }
        t
    };
    let merge = |mut s: Tables, t: Tables| -> Tables {
        for (sk, tk) in s.iter_mut().zip(t.iter()) {
            sk.iter_mut().zip(tk).for_each(|(a, b)| *a ^= b);
        }
        s
    };

    cfg_if::cfg_if! {
        if #[cfg(feature = "multithread")] {
            let t = a.par_chunks(CHUNK).zip(b.par_chunks(CHUNK)).map(partial).reduce(new_tables, merge);
        } else {
            let t = a.chunks(CHUNK).zip(b.chunks(CHUNK)).map(partial).fold(new_tables(), merge);
        }
    }

    // row 8k + l of aᵀb: ⊕ t[k][x] over the bytes x with bit l set
    from_fn(|r| {
        let (k, l) = (r / 8, r % 8);
        (0..256).filter(|x| (x >> l) & 1 == 1).fold(0, |acc, x| acc ^ t[k][x])
    })
}

/// `out ^= v·m`, for a block `v` (one `u64` per row).
fn mul_acc(out: &mut [u64], v: &[u64], m: &Mat64) {
    assert_eq!(out.len(), v.len());
    let t = byte_tables(m);
    let apply = |w: u64| t.iter().enumerate().fold(0, |acc, (k, tk)| acc ^ tk[((w >> (8 * k)) & 0xff) as usize]);

    cfg_if::cfg_if! {
        if #[cfg(feature = "multithread")] {
            out.par_iter_mut().zip(v.par_iter()).for_each(|(o, &w)| *o ^= apply(w));
        } else {
            out.iter_mut().zip(v.iter()).for_each(|(o, &w)| *o ^= apply(w));
        }
    }
}

type Tables = Box<[[u64; 256]; 8]>;

fn new_tables() -> Tables {
    Box::new([[0; 256]; 8])
}

/// `t[k][x] = (x << 8k)·m`, so that `w·m` takes 8 lookups.
fn byte_tables(m: &Mat64) -> Tables {
    let mut t = new_tables();
    for (k, tk) in t.iter_mut().enumerate() {
        for x in 1..256usize {
            tk[x] = tk[x & (x - 1)] ^ m[8 * k + x.trailing_zeros() as usize];
        }
    }
    t
}

/// A sparse 0/1 matrix in CSR form.
struct F2Mat {
    n_rows: usize,
    n_cols: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<u32>,
}

impl F2Mat {
    fn from_entries(n_rows: usize, n_cols: usize, entries: impl IntoIterator<Item = (usize, usize)>) -> Self {
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

    fn from_rows(n_cols: usize, rows: Vec<Vec<u32>>) -> Self {
        let row_ptr = once(0).chain(rows.iter().scan(0, |acc, r| {
            *acc += r.len();
            Some(*acc)
        })).collect();
        Self { n_rows: rows.len(), n_cols, row_ptr, col_idx: rows.concat() }
    }

    fn transpose(&self) -> Self {
        let entries = (0..self.n_rows).flat_map(|i| self.row(i).iter().map(move |&j| (j as usize, i)));
        Self::from_entries(self.n_cols, self.n_rows, entries)
    }

    fn row(&self, i: usize) -> &[u32] {
        &self.col_idx[self.row_ptr[i]..self.row_ptr[i + 1]]
    }

    /// The block product `self·v`, for a block `v` with one `u64` per column.
    fn mul(&self, v: &[u64]) -> Vec<u64> {
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

/// A small deterministic generator for the random start blocks.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let z = (self.0 ^ (self.0 >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        let z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use yui_matrix::sparse::pluq::solve_pluq;

    const SHAPES: [(usize, usize); 3] = [(400, 300), (300, 400), (1000, 800)];

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
