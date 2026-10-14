use std::array::from_fn;

use log::debug;

use super::f2mat::F2Mat;
use super::mat64::{Mat64, add_id, inner, mul64, mul_acc};
use super::mix::{mix_rows, unmix};
use super::rng::SplitMix64;

// Block Lanczos after [Mon95] P. L. Montgomery, "A Block Lanczos Algorithm for Finding
// Dependencies over GF(2)", EUROCRYPT '95, LNCS 921, pp. 106–120, 1995.
// Equation numbers (n) and sections §n below refer to [Mon95].
//
// Notation as in [Mon95]: `A = s = bᵀb` (symmetric, `n × n`) and `N = 64`. A block is an
// `n × N` matrix such as `V_i`, stored as one `u64` per row. A column subset `S_i` of the
// `N` columns is a bit mask, so that `·S_iS_iᵀ` (keeping the columns in `S_i`) is `& mask`.

/// One run of Block Lanczos on `s = bᵀb`. Returns 64 vectors of `ker b` (some may be
/// zero) as a block with one `u64` per coordinate, or `None` on a breakdown.
pub(super) fn kernel_block(b0: &F2Mat, seed: u64) -> Option<Vec<u64>> {
    // note: not in [Mon95]. Lanczos runs on `b = q·b0·r`, so that repeated rows or columns
    // of `b0` do not make `VᵀAV` vanish; see `mix_both`.
    let (b, bt, partner) = mix_both(b0, seed);

    // §7: for a random block `Y`, solve `A·X = A·Y` on the Krylov space of `V_0 = A·Y`.
    // If the run ends with `V_m = 0`, then `A·(X − Y) = 0`; in general, see the final step.
    let y = random_block(b.n_cols, seed);
    let (x, v_last) = iterate(&b, &bt, &y)?;

    // The final step of §7: vectors of `ker b` in the span of the columns of
    // `Z = [X − Y | V_m]`. Then `ker b0 = r·ker b` gives vectors of `ker b0`.
    let z_lo: Vec<u64> = x.iter().zip(&y).map(|(a, b)| a ^ b).collect();
    let w = unmix(&combine(&b, &z_lo, &v_last)?, &partner);

    // an exact check against b0 (it holds by construction)
    let bad = b0.mul(&w).iter().fold(0, |acc, r| acc | r);
    if bad != 0 {
        debug!("lanczos: {} vectors failed the check against the unmixed matrix", bad.count_ones());
    }
    Some(w.iter().map(|r| r & !bad).collect())
}

/// `b = q·b0·r` for random unit-triangular `q`, `r` (see `mix_rows`), with `bᵀ` and the
/// partners of `r`. Repeated rows or columns of `b0` do not cancel in `bᵀb`, and `ker b0 = r·ker b`.
fn mix_both(b0: &F2Mat, seed: u64) -> (F2Mat, F2Mat, Vec<Option<usize>>) {
    let qb0 = mix_rows(b0, seed).0;
    let (bt, partner) = mix_rows(&qb0.transpose(), !seed);
    (bt.transpose(), bt, partner)
}

fn random_block(n: usize, seed: u64) -> Vec<u64> {
    let mut rng = SplitMix64(seed);
    (0..n).map(|_| rng.next()).collect()
}

/// Runs the recurrence on `s = bᵀb` from `V_0 = s·y`, until `V_mᵀ s V_m = 0` or the column
/// selection fails. Returns `X` and `V_m`, or `None` if it does not converge.
fn iterate(b: &F2Mat, bt: &F2Mat, y: &[u64]) -> Option<(Vec<u64>, Vec<u64>)> {
    let n = b.n_cols;
    let s = |v: &[u64]| bt.mul(&b.mul(v));
    let mut rec = Recurrence::new(s(y)); // solving s·x = s·y gives x − y ∈ ker s

    // note: about n / (N − 0.76) iterations are expected (abstract, §8.1), hence n / 63 below.
    let max_iter = n / 60 + 100;
    let mut iter = 0;
    let mut dim = 0; // Σ |S_i|, the dimension of the explored Krylov space

    loop {
        if iter > max_iter {
            debug!("lanczos: no convergence after {iter} iterations");
            return None;
        }

        match rec.step(s) {
            Ok(mask) => {
                dim += mask.count_ones() as usize;
            }
            Err(Stop::Converged) => break,
            Err(Stop::NoSelection) => {
                debug!("lanczos: column selection failed at iter {iter}, stopping");
                break;
            }
        }

        iter += 1;
        if iter % 1000 == 0 {
            debug!("lanczos: iter {iter} (≤ {} expected)", n / 63);
        }
    }

    debug!("lanczos: converged after {iter} iterations, explored dimension {dim}");
    let [v_last, _, _] = rec.v;
    Some((rec.x, v_last))
}

/// The state of the recurrence after step `i`. The blocks `W_j = V_j S_j` (`j < i`) are
/// pairwise `A`-orthogonal, `W_jᵀ A W_k = 0` for `j ≠ k` ((7), (11)), and span the explored Krylov space.
struct Recurrence {
    v0: Vec<u64>,          // V_0
    v: [Vec<u64>; 3],      // V_i, V_{i-1}, V_{i-2}
    winv: [Mat64; 2],      // W_{i-1}^inv, W_{i-2}^inv
    vt_sv_prev: Mat64,     // V_{i-1}ᵀ s V_{i-1}
    vt_s2v_prev: Mat64,    // V_{i-1}ᵀ s² V_{i-1}
    mask_prev: u64,        // S_{i-1}, as a column mask
    x: Vec<u64>,           // Σ_{j < i} V_j W_j^inv V_jᵀ V_0
}

enum Stop {
    Converged,
    NoSelection,
}

impl Recurrence {
    /// The state before step 0, with `V_j = 0`, `W_j^inv = 0` and `S_j = I_N` for `j < 0`
    /// as defined after (19).
    fn new(v0: Vec<u64>) -> Self {
        let n = v0.len();
        Self {
            v: [v0.clone(), vec![0; n], vec![0; n]],
            v0,
            winv: [[0; 64]; 2],
            vt_sv_prev: [0; 64],
            vt_s2v_prev: [0; 64],
            mask_prev: u64::MAX,
            x: vec![0; n],
        }
    }

    /// Computes `V_{i+1}`, adds `V_i W_i^inv V_iᵀ V_0` to `X`, and returns `S_i`.
    fn step(&mut self, s: impl Fn(&[u64]) -> Vec<u64>) -> Result<u64, Stop> {
        // T_i = V_iᵀ A V_i and V_iᵀ A² V_i. The run ends when T_i = 0 (§7).
        let v = &self.v;
        let sv = s(&v[0]);
        let vt_sv = inner(&v[0], &sv);
        if vt_sv == [0; 64] {
            return Err(Stop::Converged);
        }
        let vt_s2v = inner(&sv, &sv);

        // S_i and W_i^inv = S_i (S_iᵀ T_i S_i)⁻¹ S_iᵀ (16), by Fig. 1 (see `select_columns`).
        // note: near the end, T_i may be too small to cover the columns skipped at the
        // previous step. We stop there; `combine` checks every vector it returns.
        let nonzero = v[0].iter().fold(0, |acc, r| acc | r); // the nonzero columns of V_i
        let Some((mask, winv0)) = select_columns(&vt_sv, self.mask_prev, nonzero) else {
            return Err(Stop::NoSelection);
        };

        // The coefficients (19), over 𝔽₂ (− = +):
        //   D_{i+1} = I − W_i^inv (V_iᵀA²V_i S_iS_iᵀ + V_iᵀAV_i),
        //   E_{i+1} = − W_{i-1}^inv V_iᵀAV_i S_iS_iᵀ,
        //   F_{i+1} = − W_{i-2}^inv (I − V_{i-1}ᵀAV_{i-1} W_{i-1}^inv)
        //               (V_{i-1}ᵀA²V_{i-1} S_{i-1}S_{i-1}ᵀ + V_{i-1}ᵀAV_{i-1}) S_iS_iᵀ.
        let (winv, mask_prev) = (&self.winv, self.mask_prev);
        let d = add_id(&mul64(&winv0, &from_fn(|i| (vt_s2v[i] & mask) ^ vt_sv[i])));
        let e = mul64(&winv[0], &vt_sv).map(|r| r & mask);
        let f = mul64(
            &mul64(&winv[1], &add_id(&mul64(&self.vt_sv_prev, &winv[0]))),
            &from_fn(|i| ((self.vt_s2v_prev[i] & mask_prev) ^ self.vt_sv_prev[i]) & mask),
        );

        // The recurrence (18): V_{i+1} = A V_i S_iS_iᵀ + V_i D_{i+1} + V_{i-1} E_{i+1} + V_{i-2} F_{i+1}.
        let mut v_next: Vec<u64> = sv.iter().map(|r| r & mask).collect();
        mul_acc(&mut v_next, &v[0], &d);
        mul_acc(&mut v_next, &v[1], &e);
        mul_acc(&mut v_next, &v[2], &f);

        // X = Σ_i V_i W_i^inv V_iᵀ V_0 (20), which solves A·X = V_0 on the span of the W_i (8).
        mul_acc(&mut self.x, &v[0], &mul64(&winv0, &inner(&v[0], &self.v0)));

        // shift: i → i + 1
        self.v.rotate_right(1);
        self.v[0] = v_next;
        self.winv = [winv0, self.winv[0]];
        self.vt_sv_prev = vt_sv;
        self.vt_s2v_prev = vt_s2v;
        self.mask_prev = mask;
        Ok(mask)
    }
}

/// Montgomery's choice of `S_i` (§8, Fig. 1): a column mask `S` and `W^inv = S (Sᵀ t S)⁻¹ Sᵀ` (16),
/// where `S` must contain every column not in `S_{i-1}` (15) among the `nonzero` columns of `V_i`.
/// `None` on failure (a breakdown).
fn select_columns(t: &Mat64, mask_prev: u64, nonzero: u64) -> Option<(u64, Mat64)> {
    // Fig. 1: Gaussian elimination on [T | I], with the columns in S_{i-1} last,
    // so that the others are chosen for S_i whenever T allows.
    let mut m: [(u64, u64); 64] = from_fn(|i| (t[i], 1 << i)); // [t | I]
    let order: Vec<usize> = (0..64).filter(|&c| (mask_prev >> c) & 1 == 0)
        .chain((0..64).filter(|&c| (mask_prev >> c) & 1 == 1))
        .collect();
    let mut mask = 0;

    for i in 0..64 {
        let c = order[i];
        let bit = 1 << c;

        if let Some(j) = (i..64).find(|&j| m[order[j]].0 & bit != 0) {
            // a pivot in the T half: c ∈ S
            m.swap(c, order[j]);
            eliminate(&mut m, c, |row| row.0 & bit != 0);
            mask |= bit;
        } else {
            // no pivot in the T half: c ∉ S. Eliminate with the I half, and drop row c.
            let Some(j) = (i..64).find(|&j| m[order[j]].1 & bit != 0) else {
                debug!("lanczos: select_columns: no pivot in the inverse half");
                return None;
            };
            m.swap(c, order[j]);
            eliminate(&mut m, c, |row| row.1 & bit != 0);
            m[c] = (0, 0);
        }
    }

    // Every nonzero column of V_i must be in S_i ∪ S_{i-1}, by (15), as checked in §8.
    let uncovered = !(mask | mask_prev) & nonzero;
    if mask == 0 || uncovered != 0 {
        debug!("lanczos: select_columns: |S| = {}, |S_prev| = {}, uncovered = {}",
            mask.count_ones(), mask_prev.count_ones(), uncovered.count_ones());
        return None;
    }

    // the I half now holds W^inv
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

/// Montgomery's final step (§7): with `Z = [X − Y | V_m]`, returns `Z·u` for the `u` with `b·Z·u = 0`.
fn combine(b: &F2Mat, z_lo: &[u64], z_hi: &[u64]) -> Option<Vec<u64>> {
    // b·Z, an m × 2N matrix (one u128 per row), and up to N of its null vectors u ∈ 𝔽₂^{2N}.
    // note: §7 keeps all of them (up to 2N); we keep N, one block.
    let (bz_lo, bz_hi) = (b.mul(z_lo), b.mul(z_hi));
    let mut rows: Vec<u128> = bz_lo.iter().zip(&bz_hi).map(|(&l, &h)| l as u128 | (h as u128) << 64).collect();
    let tags: Vec<u128> = null_space_128(&mut rows).into_iter().take(64).collect();

    if tags.is_empty() {
        return None;
    }

    // W = Z·U = (X − Y)·U_lo + V_m·U_hi, where U is the 2N × N matrix whose columns are
    // the null vectors, and U_lo, U_hi are its N × N halves.
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
    // tags[c]: column c of the current matrix as a combination of the original columns.
    // A column that is zero in every row at the end is a null vector, given by its tag.
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

#[cfg(test)]
mod tests {
    use std::cmp::Reverse;

    use super::*;
    use super::super::test_utils::{bit, rand_block, rand_f2mat};

    const SIZE: (usize, usize) = (600, 400);

    /// `b` with column `n - d + t` = column `2t` + column `2t + 1`, so that
    /// `e_{2t} + e_{2t+1} + e_{n-d+t}` (`t < d`) span `ker b` (when the rest has full rank).
    fn planted(rng: &mut SplitMix64, (m, n): (usize, usize), d: usize) -> F2Mat {
        let base = rand_f2mat(rng, (m, n - d), 0.02);
        let mut entries: Vec<(usize, usize)> = (0..m).flat_map(|i| {
            base.row(i).iter().flat_map(move |&j| {
                let j = j as usize;
                [Some((i, j)), (j < 2 * d).then_some((i, n - d + j / 2))]
            }).flatten()
        }).collect();
        entries.sort_unstable();

        // note: (i, n-d+t) appears twice when both columns 2t and 2t+1 hit row i; they cancel.
        let entries = entries.chunk_by(|a, b| a == b).filter(|c| c.len() % 2 == 1).map(|c| c[0]);
        F2Mat::from_entries(m, n, entries.collect::<Vec<_>>())
    }

    /// The column rank of a block, i.e. the rank of its rows as 64-bit vectors.
    fn rank_64(w: &[u64]) -> usize {
        rank_128(w.iter().map(|&r| r as u128).collect())
    }

    /// The rank over `𝔽₂` of the rows, kept with distinct leading bits in decreasing order.
    fn rank_128(rows: Vec<u128>) -> usize {
        let lead = |x: u128| 127 - x.leading_zeros();
        rows.into_iter().fold(Vec::<u128>::new(), |mut pivots, r| {
            let x = pivots.iter().fold(r, |x, &p| if (x >> lead(p)) & 1 == 1 { x ^ p } else { x });
            if x != 0 {
                pivots.push(x);
                pivots.sort_by_key(|&p| Reverse(lead(p)));
            }
            pivots
        }).len()
    }

    #[test]
    fn null_space_128_is_the_kernel() {
        let mut rng = SplitMix64(17);
        // columns 100..128 are zero, and column 99 = column 0 + column 1.
        let low = (1u128 << 99) - 1;
        for n_rows in [0, 50, 127, 200] {
            let rows0: Vec<u128> = (0..n_rows).map(|_| {
                let r = (rng.next() as u128 | (rng.next() as u128) << 64) & low;
                r | ((r ^ (r >> 1)) & 1) << 99
            }).collect();

            let tags = null_space_128(&mut rows0.clone());
            assert!(tags.iter().all(|&t| rows0.iter().all(|&r| (r & t).count_ones() % 2 == 0)), "null vectors");
            assert_eq!(rank_128(tags.clone()), tags.len(), "independent");
            assert_eq!(tags.len(), 128 - rank_128(rows0), "a basis");
        }
    }

    #[test]
    fn select_columns_inverts_on_the_selection() {
        let mut rng = SplitMix64(18);
        let id: Mat64 = from_fn(|i| 1 << i);

        // t invertible: S is everything, and W^inv = t⁻¹.
        let v = rand_block(&mut rng, 2000);
        let t = inner(&v, &v);
        let (mask, winv) = select_columns(&t, 0, u64::MAX).expect("t is invertible");
        assert_eq!(mask, u64::MAX);
        assert_eq!(mul64(&winv, &t), id);

        // t supported on columns 0..40 and S_prev = 40..64: S = 0..40.
        let low = (1u64 << 40) - 1;
        let v: Vec<u64> = rand_block(&mut rng, 2000).iter().map(|x| x & low).collect();
        let t = inner(&v, &v);
        let (mask, winv) = select_columns(&t, !low, u64::MAX).expect("a selection exists");
        let wt = mul64(&winv, &t);
        assert_eq!(mask, low);
        assert_eq!(mul64(&wt, &winv), winv);
        assert!((0..64).all(|i| winv[i] & !mask == 0 && (bit(mask, i) || winv[i] == 0)));
        assert!((0..40).all(|i| wt[i] == 1 << i));
    }

    // §8: an uncovered column outside S_{i-1} is a breakdown only if it is nonzero in V_i.
    #[test]
    fn select_columns_ignores_zero_columns() {
        // t = I on columns 0..30 and S_prev = 40..64: columns 30..40 cannot be selected.
        // note: v has the unit rows, so that t is invertible there (a random vᵀv often is not).
        let (low, sel) = ((1u64 << 40) - 1, (1u64 << 30) - 1);
        let v: Vec<u64> = (0..30).map(|k| 1 << k).collect();
        let t = inner(&v, &v);

        let (mask, _) = select_columns(&t, !low, sel).expect("columns 30..40 are zero");
        assert_eq!(mask, sel);
        assert!(select_columns(&t, !low, u64::MAX).is_none(), "columns 30..40 are nonzero");
    }

    // Montgomery's invariant: the blocks W_i = V_i S_i are pairwise s-orthogonal.
    #[test]
    fn recurrence_blocks_are_s_orthogonal() {
        let mut rng = SplitMix64(19);
        let b = rand_f2mat(&mut rng, SIZE, 0.02);
        let bt = b.transpose();
        let s = |v: &[u64]| bt.mul(&b.mul(v));
        let mut rec = Recurrence::new(s(&rand_block(&mut rng, b.n_cols)));

        let mut ws: Vec<Vec<u64>> = vec![];
        loop {
            let v = rec.v[0].clone();
            let Ok(mask) = rec.step(s) else {
                break;
            };
            ws.push(v.iter().map(|r| r & mask).collect());
        }

        assert!(ws.len() > 2, "too few steps: {}", ws.len());
        for (i, wi) in ws.iter().enumerate() {
            let swi = s(wi);
            assert!(ws[..i].iter().all(|wj| inner(wj, &swi) == [0; 64]), "W_j s W_{i} ≠ 0");
        }
    }

    #[test]
    fn iterate_and_combine_find_kernel_vectors() {
        let mut rng = SplitMix64(20);
        for d in [1, 10, 40] {
            let b = planted(&mut rng, SIZE, d);
            let y = rand_block(&mut rng, b.n_cols);
            let (x, v_last) = iterate(&b, &b.transpose(), &y).expect("converges");

            let z_lo: Vec<u64> = x.iter().zip(&y).map(|(a, b)| a ^ b).collect();
            let w = combine(&b, &z_lo, &v_last).expect("kernel vectors");
            assert!(b.mul(&w).iter().all(|&r| r == 0), "d = {d}: not in ker b");
            assert!((1..=d).contains(&rank_64(&w)), "d = {d}: rank {}", rank_64(&w));
        }
    }

    // With Z = [z_lo | z_hi], combine returns min(128 − rank(b·Z), 64) independent vectors of ker b.
    #[test]
    fn combine_returns_kernel_vectors_in_the_span() {
        let mut rng = SplitMix64(21);
        let (d, n) = (10, SIZE.1);
        let b = planted(&mut rng, SIZE, d);

        // columns 0..10 of z_hi: random combinations of the planted kernel basis
        let z_lo = rand_block(&mut rng, n);
        let mut z_hi = rand_block(&mut rng, n).iter().map(|r| r & !0x3ff).collect::<Vec<u64>>();
        for c in 0..10 {
            for t in (0..d).filter(|_| rng.next() & 1 == 1) {
                [2 * t, 2 * t + 1, n - d + t].iter().for_each(|&r| z_hi[r] ^= 1 << c);
            }
        }

        let bz: Vec<u128> = b.mul(&z_lo).iter().zip(b.mul(&z_hi))
            .map(|(&l, h)| l as u128 | (h as u128) << 64)
            .collect();
        let expected = (128 - rank_128(bz)).min(64);

        let w = combine(&b, &z_lo, &z_hi).expect("kernel vectors");
        assert!(b.mul(&w).iter().all(|&r| r == 0));
        assert_eq!(rank_64(&w), expected);
        assert!(expected >= 1);
    }
}
