use ahash::AHashMap;
use crate::grid_diagram::{GridDiagram, GridRect};
use crate::grid_complex::GridComplexGen;

pub type GeneratorTable = AHashMap<(isize, isize), Vec<GridComplexGen>>;

pub fn produce(diagram: &GridDiagram) -> GeneratorTable {
    let n = diagram.grid_number();
    let size = 2 * n;
    let trans = heap_transpositions(n);

    let mut result: GeneratorTable = AHashMap::new();
    let mut seq: Vec<u8> = (0..n).collect();
    let mut x = GridComplexGen::new(diagram, &seq);
    let mut pts = x.points();

    insert(&mut result, x.clone());

    for &(i, j) in trans.iter() {
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

        insert(&mut result, x.clone());
    }

    result
}

fn insert(table: &mut GeneratorTable, x: GridComplexGen) {
    table.entry((x.maslov_deg, x.alexander_deg)).or_default().push(x);
}

// Heap's algorithm: https://en.wikipedia.org/wiki/Heap%27s_algorithm
fn heap_transpositions(n: u8) -> Vec<(u8, u8)> {
    let capacity = (1..=(n as usize)).product::<usize>(); // n!
    let mut result = Vec::with_capacity(capacity);

    fn generate(k: u8, result: &mut Vec<(u8, u8)>) {
        if k <= 1 { return; }
        generate(k - 1, result);
        for i in 0..k - 1 {
            let t = if k % 2 == 0 { (i, k - 1) } else { (0, k - 1) };
            result.push(t);
            generate(k - 1, result);
        }
    }

    generate(n, &mut result);
    result
}

use crate::grid_diagram::GridPoint;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heap_transpositions_count() {
        // n! - 1 transpositions for n elements
        assert_eq!(heap_transpositions(3).len(), 5);  // 3! - 1
        assert_eq!(heap_transpositions(4).len(), 23); // 4! - 1
    }

    #[test]
    fn test_produce_trefoil() {
        let diagram = GridDiagram::from_sequences(&[1, 2, 0], &[0, 1, 2]);
        let n = diagram.grid_number();
        let table = produce(&diagram);

        // 3x3 grid has 3! = 6 generators
        let total: usize = table.values().map(|v| v.len()).sum();
        assert_eq!(total, 6);

        // Check all generators are produced correctly
        let trans = heap_transpositions(3);
        let mut seq: Vec<u8> = (0..n).collect();

        let x = GridComplexGen::new(&diagram, &seq);
        assert!(table[&(x.maslov_deg, x.alexander_deg)].contains(&x));

        for &(i, j) in trans.iter() { 
            seq.swap(i as usize, j as usize);
            let x = GridComplexGen::new(&diagram, &seq);
            assert!(table[&(x.maslov_deg, x.alexander_deg)].contains(&x));
        }
    }
}
