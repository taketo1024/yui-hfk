// Heap's algorithm (iterative): https://en.wikipedia.org/wiki/Heap%27s_algorithm
pub struct HeapTranspositions {
    n: u8,
    c: Vec<u8>,
    i: u8,
}

impl HeapTranspositions {
    pub fn new(n: u8) -> Self {
        HeapTranspositions { n, c: vec![0; n as usize], i: 1 }
    }
}

impl Iterator for HeapTranspositions {
    type Item = (u8, u8);

    fn next(&mut self) -> Option<(u8, u8)> {
        while self.i < self.n {
            if self.c[self.i as usize] < self.i {
                let swap = if self.i % 2 == 0 {
                    (0, self.i)
                } else {
                    (self.c[self.i as usize], self.i)
                };
                self.c[self.i as usize] += 1;
                self.i = 1;
                return Some(swap);
            } else {
                self.c[self.i as usize] = 0;
                self.i += 1;
            }
        }
        None
    }
}

pub fn heap_transpositions(n: u8) -> HeapTranspositions {
    HeapTranspositions::new(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_heap_transpositions_entries() {
        let v2: Vec<_> = heap_transpositions(2).collect();
        assert_eq!(v2, vec![(0, 1)]);
        let v3: Vec<_> = heap_transpositions(3).collect();
        assert_eq!(v3, vec![(0, 1), (0, 2), (0, 1), (0, 2), (0, 1)]);
    }

    #[test]
    fn test_heap_transpositions_generates_all_permutations() {
        // Applying the transpositions to the identity must visit every
        // permutation of {0..n-1} exactly once (including the start).
        for n in 2u8..=5 {
            let trans = heap_transpositions(n);
            let mut seq: Vec<u8> = (0..n).collect();
            let mut seen: HashSet<Vec<u8>> = HashSet::new();

            seen.insert(seq.clone());

            for (i, j) in trans {
                seq.swap(i as usize, j as usize);
                assert!(seen.insert(seq.clone()), "duplicate permutation at n={n}: {seq:?}");
            }
            let n_fact: usize = (1..=n as usize).product();
            assert_eq!(seen.len(), n_fact, "wrong count at n={n}");
        }
    }
}
