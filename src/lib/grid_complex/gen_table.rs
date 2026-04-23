use ahash::AHashSet;
use crate::grid_diagram::{GridDiagram, GridPoint, GridRect};
use crate::misc::heap::heap_transpositions;
use super::gen::GridComplexGen;

pub struct GridComplexGenTable {
    inner: AHashSet<GridComplexGen>,
}

impl GridComplexGenTable {
    fn new() -> Self {
        GridComplexGenTable { inner: AHashSet::new() }
    }

    pub fn generate_from(diagram: &GridDiagram) -> Self { 
        Self::generate_filtered(diagram, |_| true)
    }

    pub fn generate_filtered(diagram: &GridDiagram, pred: impl Fn(&GridComplexGen) -> bool) -> Self { 
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

            if pred(&x) {
                result.insert(x.clone());
            }
        }

        result
    }

    pub fn len(&self) -> usize { 
        self.inner.len()
    }

    pub fn insert(&mut self, x: GridComplexGen) {
        self.inner.insert(x);
    }

    pub fn contains(&self, x: &GridComplexGen) -> bool {
        self.inner.contains(x)
    }

    pub fn iter(&self) -> impl Iterator<Item = &GridComplexGen> {
        self.inner.iter()
    }

    pub fn gen_from_seq(&self, seq: &[u8]) -> Option<&GridComplexGen> {
        let code = GridComplexGen::encode(seq);
        self.inner.iter().find(|g| g.code == code)
    }

    pub fn gens_of_degree(&self, m: isize, a: isize) -> impl Iterator<Item = &GridComplexGen> {
        self.inner.iter().filter(move |g| g.maslov_deg == m && g.alexander_deg == a)
    }

    pub fn support(&self) -> impl Iterator<Item = (isize, isize)> {
        let mut degs: AHashSet<(isize, isize)> = AHashSet::new();
        for g in &self.inner {
            degs.insert((g.maslov_deg, g.alexander_deg));
        }
        degs.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_produce() {
        let diagram = GridDiagram::from_sequences(&[1, 2, 0], &[0, 1, 2]);
        let n = diagram.grid_number();
        let table = GridComplexGenTable::generate_from(&diagram);

        // 3x3 grid has 3! = 6 generators
        assert_eq!(table.iter().count(), 6);

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

    #[test]
    fn test_produce_trefoil() {
        let diagram = GridDiagram::trefoil();
        let table = GridComplexGenTable::generate_from(&diagram);

        let n = diagram.grid_number();
        assert_eq!(n, 5);
        assert_eq!(table.iter().count(), 120);

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
