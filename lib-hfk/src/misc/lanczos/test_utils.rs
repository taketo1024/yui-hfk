use super::f2mat::F2Mat;
use super::rng::SplitMix64;

pub(super) const SHAPES: [(usize, usize); 3] = [(400, 300), (300, 400), (1000, 800)];

pub(super) fn rand_f2mat(rng: &mut SplitMix64, (m, n): (usize, usize), density: f64) -> F2Mat {
    let threshold = (density * u64::MAX as f64) as u64;
    let entries: Vec<(usize, usize)> = (0..m).flat_map(|i| (0..n).map(move |j| (i, j)))
        .filter(|_| rng.next() < threshold)
        .collect();
    F2Mat::from_entries(m, n, entries)
}

pub(super) fn dense(a: &F2Mat) -> Vec<Vec<bool>> {
    (0..a.n_rows).map(|i| {
        let mut row = vec![false; a.n_cols];
        a.row(i).iter().for_each(|&j| row[j as usize] ^= true);
        row
    }).collect()
}

pub(super) fn rand_block(rng: &mut SplitMix64, n: usize) -> Vec<u64> {
    (0..n).map(|_| rng.next()).collect()
}

pub(super) fn bit(x: u64, i: usize) -> bool {
    (x >> i) & 1 == 1
}
