use std::array::from_fn;

#[cfg(feature = "multithread")]
use rayon::prelude::*;

/// A 64 × 64 matrix over `𝔽₂`: row `i` is a word, bit `j` is column `j`.
pub(super) type Mat64 = [u64; 64];

pub(super) fn mul64(a: &Mat64, b: &Mat64) -> Mat64 {
    from_fn(|i| (0..64).filter(|&j| (a[i] >> j) & 1 == 1).fold(0, |acc, j| acc ^ b[j]))
}

pub(super) fn add_id(a: &Mat64) -> Mat64 {
    from_fn(|i| a[i] ^ (1 << i))
}

/// For blocks `a`, `b` (one `u64` per row), returns `aᵀb`.
pub(super) fn inner(a: &[u64], b: &[u64]) -> Mat64 {
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
pub(super) fn mul_acc(out: &mut [u64], v: &[u64], m: &Mat64) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::rng::SplitMix64;
    use super::super::test_utils::{bit, rand_block};

    #[test]
    fn mul64_and_add_id() {
        let mut rng = SplitMix64(13);
        let id: Mat64 = from_fn(|i| 1 << i);
        for _ in 0..10 {
            let a: Mat64 = from_fn(|_| rng.next());
            let b: Mat64 = from_fn(|_| rng.next());
            let ab = mul64(&a, &b);

            assert!((0..64).all(|i| (0..64).all(|j| {
                bit(ab[i], j) == (0..64).fold(false, |acc, k| acc ^ (bit(a[i], k) && bit(b[k], j)))
            })));
            assert_eq!(mul64(&a, &id), a);
            assert_eq!(mul64(&id, &a), a);
            assert!((0..64).all(|i| add_id(&a)[i] ^ a[i] == 1 << i));
        }
    }

    #[test]
    fn inner_matches_naive() {
        let mut rng = SplitMix64(14);
        // note: 40_000 is longer than the chunk size 1 << 14, so partial tables are merged.
        for len in [0, 1, 100, 40_000] {
            let (a, b) = (rand_block(&mut rng, len), rand_block(&mut rng, len));
            let expected: Mat64 = from_fn(|r| {
                a.iter().zip(&b).filter(|(x, _)| bit(**x, r)).fold(0, |acc, (_, y)| acc ^ y)
            });
            assert_eq!(inner(&a, &b), expected);
        }
    }

    #[test]
    fn mul_acc_matches_naive() {
        let mut rng = SplitMix64(15);
        let m: Mat64 = from_fn(|_| rng.next());
        let v = rand_block(&mut rng, 1000);
        let out0 = rand_block(&mut rng, 1000);

        let mut out = out0.clone();
        mul_acc(&mut out, &v, &m);

        let expected: Vec<u64> = out0.iter().zip(&v).map(|(&o, &w)| {
            (0..64).filter(|&j| bit(w, j)).fold(o, |acc, j| acc ^ m[j])
        }).collect();
        assert_eq!(out, expected);
    }

    #[test]
    fn byte_tables_match_definition() {
        let mut rng = SplitMix64(16);
        let m: Mat64 = from_fn(|_| rng.next());
        let t = byte_tables(&m);
        assert!((0..8).all(|k| (0..256usize).all(|x| {
            t[k][x] == (0..8).filter(|&l| (x >> l) & 1 == 1).fold(0, |acc, l| acc ^ m[8 * k + l])
        })));
    }
}
