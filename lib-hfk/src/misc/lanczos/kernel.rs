use std::array::from_fn;

use log::debug;

use super::f2mat::F2Mat;
use super::mat64::{Mat64, add_id, inner, mul64, mul_acc};
use super::mix::{mix_rows, unmix};
use super::rng::SplitMix64;

/// One run of Block Lanczos on `s = bᵀb`. Returns 64 vectors of `ker b` (some may be
/// zero) as a block with one `u64` per coordinate, or `None` on a breakdown.
pub(super) fn kernel_block(b0: &F2Mat, seed: u64) -> Option<Vec<u64>> {
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

#[cfg(test)]
mod tests {
    use std::cmp::Reverse;

    use super::*;
    use super::super::test_utils::{bit, rand_block};

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
        let (mask, winv) = select_columns(&t, 0).expect("t is invertible");
        assert_eq!(mask, u64::MAX);
        assert_eq!(mul64(&winv, &t), id);

        // t supported on columns 0..40 and S_prev = 40..64: S = 0..40.
        let low = (1u64 << 40) - 1;
        let v: Vec<u64> = rand_block(&mut rng, 2000).iter().map(|x| x & low).collect();
        let t = inner(&v, &v);
        let (mask, winv) = select_columns(&t, !low).expect("a selection exists");
        let wt = mul64(&winv, &t);
        assert_eq!(mask, low);
        assert_eq!(mul64(&wt, &winv), winv);
        assert!((0..64).all(|i| winv[i] & !mask == 0 && (bit(mask, i) || winv[i] == 0)));
        assert!((0..40).all(|i| wt[i] == 1 << i));
    }
}
