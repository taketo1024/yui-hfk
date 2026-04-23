use ahash::AHashMap;
use crate::grid_diagram::{GridDiagram, GridPoint, GridRect};
use crate::misc::heap::heap_transpositions;
use super::gen::GridComplexGen;

pub struct GridComplexGenTable {
    inner: AHashMap<(isize, isize), Vec<GridComplexGen>>,
}

impl GridComplexGenTable {
    fn new() -> Self {
        GridComplexGenTable { inner: AHashMap::new() }
    }

    pub fn insert(&mut self, x: GridComplexGen) {
        self.inner.entry((x.maslov_deg, x.alexander_deg)).or_default().push(x);
    }

    pub fn contains(&self, x: &GridComplexGen) -> bool {
        self.inner.get(&(x.maslov_deg, x.alexander_deg))
            .is_some_and(|v| v.contains(x))
    }

    pub fn gens_of_degree(&self, m: isize, a: isize) -> &[GridComplexGen] {
        self.inner.get(&(m, a)).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn degrees(&self) -> impl Iterator<Item = (isize, isize)> + '_ {
        self.inner.keys().copied()
    }

    pub fn all_gens(&self) -> impl Iterator<Item = &GridComplexGen> {
        self.inner.values().flatten()
    }

    fn sort_buckets(&mut self) {
        for gens in self.inner.values_mut() {
            gens.sort_unstable();
        }
    }
}

impl From<&GridDiagram> for GridComplexGenTable {
    fn from(diagram: &GridDiagram) -> Self {
        let n = diagram.grid_number();
        let size = 2 * n;
        let trans = heap_transpositions(n);

        let mut result = GridComplexGenTable::new();
        let mut seq: Vec<u8> = (0..n).collect();
        let mut x = GridComplexGen::new(diagram, &seq);
        let mut pts = x.points();

        result.insert(x.clone());

        for (i, j) in trans {
            let (i, j) = (i as usize, j as usize);
            let rect = GridRect::between(pts[i], pts[j], size);

            let n_o = rect.count_intersections(&diagram.os, size, false) as isize;
            let n_x = rect.count_intersections(&diagram.xs, size, false) as isize;
            let c   = rect.count_intersections(&pts,        size, true)  as isize;

            seq.swap(i, j);

            pts[i] = GridPoint::new(2 * i as u8, 2 * seq[i]);
            pts[j] = GridPoint::new(2 * j as u8, 2 * seq[j]);

            // See Book, p.68
            // (4.2)  M(y) - M(x) = 2 #(r ∩ Os) - 2 #(x ∩ Int(r)) - 1
            // (4.4)  A(y) - A(x) = #(r ∩ Os) - #(r ∩ Xs)

            let m_deg = x.maslov_deg + 2 * (n_o - c) - 1;
            let a_deg = x.alexander_deg + n_o - n_x;

            x = GridComplexGen::from_raw(
                GridComplexGen::encode(&seq),
                n,
                m_deg,
                a_deg
            );

            result.insert(x.clone());
        }

        result.sort_buckets();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_produce_trefoil() {
        let diagram = GridDiagram::from_sequences(&[1, 2, 0], &[0, 1, 2]);
        let n = diagram.grid_number();
        let table = GridComplexGenTable::from(&diagram);

        // 3x3 grid has 3! = 6 generators
        assert_eq!(table.all_gens().count(), 6);

        // Check all generators are produced correctly
        let mut seq: Vec<u8> = (0..n).collect();

        let x = GridComplexGen::new(&diagram, &seq);
        assert!(table.contains(&x));

        for (i, j) in heap_transpositions(n) {
            seq.swap(i as usize, j as usize);
            let x = GridComplexGen::new(&diagram, &seq);
            assert!(table.contains(&x));
        }
    }
}
