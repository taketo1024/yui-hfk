use std::ops::{Range, RangeInclusive};

use ahash::AHashMap;
use yui_homology::range_of;
use crate::grid_diagram::{GridDiagram, GridPoint, GridRect};
use crate::misc::heap::heap_transpositions;
use super::gen::GCGen;

pub struct GCGenSet {
    list: Vec<GCGen>, // sorted by (alex_deg, maslov_deg)
    dict: AHashMap<(isize, isize), Range<usize>>, // (maslov_deg, alex_deg) → range into list
    code_index: AHashMap<usize, usize>,            // code → index into list
}

impl GCGenSet {
    pub fn generate_from(diagram: &GridDiagram) -> Self {
        Self::generate_filtered(diagram, |_| true)
    }

    pub fn generate_filtered(diagram: &GridDiagram, pred: impl Fn(&GCGen) -> bool) -> Self {
        let n = diagram.grid_number();
        let size = 2 * n;
        let trans = heap_transpositions(n);

        let mut gens: Vec<GCGen> = Vec::new();
        let mut seq: Vec<u8> = (0..n).collect();
        let mut x = GCGen::new(diagram, &seq);
        let mut pts = x.points();

        if pred(&x) {
            gens.push(x.clone());
        }

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
            let a_deg = x.alex_deg + n_o - n_x;

            x = GCGen::from_raw(0, n, m_deg, a_deg);

            if pred(&x) {
                x.code = GCGen::encode(&seq);
                gens.push(x.clone());
            }
        }

        gens.sort();
        let by_deg = Self::build_index(&gens);
        let code_index = gens.iter().enumerate().map(|(i, g)| (g.code, i)).collect();

        GCGenSet { list: gens, dict: by_deg, code_index }
    }

    fn build_index(gens: &[GCGen]) -> AHashMap<(isize, isize), Range<usize>> {
        let mut by_deg = AHashMap::new();
        let mut i = 0;
        while i < gens.len() {
            let key = (gens[i].maslov_deg, gens[i].alex_deg);
            let j = i + gens[i..].partition_point(|g| g.maslov_deg == key.0 && g.alex_deg == key.1);
            by_deg.insert(key, i..j);
            i = j;
        }
        by_deg
    }

    pub fn len(&self) -> usize {
        self.list.len()
    }

    pub fn contains(&self, x: &GCGen) -> bool {
        self.raw_gens(x.maslov_deg, x.alex_deg).iter().any(|g| g.code == x.code)
    }

    pub fn iter(&self) -> impl Iterator<Item = &GCGen> {
        self.list.iter()
    }

    pub fn raw_gens(&self, maslov: isize, alex: isize) -> &[GCGen] {
        self.dict.get(&(maslov, alex)).map_or(&[], |r| &self.list[r.clone()])
    }

    pub fn gen_from_code(&self, code: usize) -> Option<&GCGen> {
        self.code_index.get(&code).map(|&i| &self.list[i])
    }

    pub fn maslov_deg_range(&self) -> RangeInclusive<isize> {
        range_of(self.dict.keys().map(|&(m, _)| m))
    }

    pub fn alex_deg_range(&self) -> RangeInclusive<isize> {
        range_of(self.dict.keys().map(|&(_, a)| a))
    }

    pub fn has_gen_at(&self, maslov: isize, alex: isize) -> bool {
        self.dict.contains_key(&(maslov, alex))
    }

    pub fn support(&self) -> impl Iterator<Item = &(isize, isize)> { 
        self.dict.keys()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_produce() {
        let diagram = GridDiagram::from_sequences(&[1, 2, 0], &[0, 1, 2]);
        let n = diagram.grid_number();
        let table = GCGenSet::generate_from(&diagram);

        // 3x3 grid has 3! = 6 generators
        assert_eq!(table.iter().count(), 6);

        // Check all generators are produced correctly
        let mut seq: Vec<u8> = (0..n).collect();

        let x = GCGen::new(&diagram, &seq);
        assert!(table.contains(&x));

        for (i, j) in heap_transpositions(n) {
            seq.swap(i as usize, j as usize);
            let x = GCGen::new(&diagram, &seq);
            assert!(table.contains(&x));
        }
    }

    #[test]
    fn test_produce_trefoil() {
        let diagram = GridDiagram::trefoil();
        let table = GCGenSet::generate_from(&diagram);

        let n = diagram.grid_number();
        assert_eq!(n, 5);
        assert_eq!(table.iter().count(), 120);

        // Check all generators are produced correctly
        let mut seq: Vec<u8> = (0..n).collect();

        let x = GCGen::new(&diagram, &seq);
        assert!(table.contains(&x));

        for (i, j) in heap_transpositions(n) {
            seq.swap(i as usize, j as usize);
            let x = GCGen::new(&diagram, &seq);
            assert!(table.contains(&x));
        }
    }

    #[test]
    fn test_raw_gens() {
        let diagram = GridDiagram::trefoil();
        let table = GCGenSet::generate_from(&diagram);

        for (&(m, a), range) in &table.dict {
            assert!(!range.is_empty());
            for g in &table.list[range.clone()] {
                assert_eq!(g.maslov_deg, m);
                assert_eq!(g.alex_deg, a);
            }
        }

        // Total across all bidegrees equals total count
        let total: usize = table.dict.values().map(|r| r.len()).sum();
        assert_eq!(total, table.len());
    }
}
