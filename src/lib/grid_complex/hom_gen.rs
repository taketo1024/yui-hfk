use num_traits::One;
use yui_core::{lc::Lc, num::FF2, poly::MultiVar};

use crate::grid_diagram::{GridPoint, GridRect};
use super::{GCFiltered, GCGen, GCInfGen, GridComplex};

impl GridComplex<GCFiltered> {
    /// Computes the homological generator z for the filtered grid complex.
    ///
    /// Starting from the top-left generator, applies a sequence of pentagon maps
    /// (one per adjacent transposition in the bubble sort from O₀ to O₁) to
    /// produce a cycle in the filtered complex.
    pub fn hom_gen(&self) -> Lc<GCInfGen, FF2> {
        let n = self.diagram.grid_number() as usize;
        let grid_size = self.diagram.size;

        // curr_os[i] = row index of the O-mark in column i
        let mut curr_os: Vec<u8> = self.diagram.os.iter().map(|p| (p.y - 1) / 2).collect();

        // O0 = [n-1, …, 0] (descending); O1 = actual O-row arrangement
        let o0: Vec<u8> = (0..n as u8).rev().collect();
        let o1 = curr_os.clone();

        let trans = adjacent_transpositions(&o0, &o1);

        // TL generator: seq[i] = (O0[i] + 1) % n = (n − i) % n
        let tl_seq: Vec<u8> = o0.iter().map(|&o| ((o as usize + 1) % n) as u8).collect();
        let x0 = GCGen::from_raw(GCGen::encode(&tl_seq), n as u8, 0, 0);

        let mut z: Lc<GCInfGen, FF2> = Lc::from(GCInfGen::new(x0, MultiVar::<'U', usize>::one()));

        curr_os = o0;
        for (i, j) in trans {
            z = apply_pentagon(z, i, j, &curr_os, grid_size);
            curr_os.swap(i, j);
        }

        // Replace dummy generators with actual generators from the diagram
        z.into_iter().filter_map(|(inf_gen, coeff)| {
            let seq = inf_gen.gen.sequence();
            self.gens.gen_from_seq(&seq).map(|actual| {
                (GCInfGen::new(actual.clone(), inf_gen.u), coeff)
            })
        }).collect()
    }
}

/// Bubble-sort `from` into `to`, recording each adjacent swap (k, k+1).
fn adjacent_transpositions(from: &[u8], to: &[u8]) -> Vec<(usize, usize)> {
    let n = from.len();
    let mut curr = from.to_vec();
    let mut trans = Vec::new();
    for target in 0..n {
        let j = curr.iter().position(|&x| x == to[target]).unwrap();
        for k in (target..j).rev() {
            curr.swap(k, k + 1);
            trans.push((k, k + 1));
        }
    }
    trans
}

/// Apply the pentagon map for transposition (i, j) to z.
fn apply_pentagon(
    z: Lc<GCInfGen, FF2>,
    i: usize,
    j: usize,
    curr_os: &[u8],
    grid_size: u8,
) -> Lc<GCInfGen, FF2> {
    let n = curr_os.len();

    // Oi, Oj: odd-coordinate O-mark positions; v: pentagon vertex (one left of Oj)
    let oi = GridPoint::new(2 * i as u8 + 1, 2 * curr_os[i] + 1);
    let oj = GridPoint::new(2 * j as u8 + 1, 2 * curr_os[j] + 1);
    let v  = GridPoint::new(2 * j as u8, oj.y);

    let o_pts: Vec<GridPoint> = (0..n)
        .map(|r| GridPoint::new(2 * r as u8 + 1, 2 * curr_os[r] + 1))
        .collect();

    z.into_iter()
        .flat_map(|(inf_gen, _)| {
            let x_seq = inf_gen.gen.sequence();
            let x_pts = inf_gen.gen.points();
            let q = swap_u(&inf_gen.u, i, j, n);
            let o_ref = o_pts.as_slice();
            (0..n).filter_map(move |k| {
                if k == j { return None; }
                let mut y_seq = x_seq.clone();
                y_seq.swap(j, k);
                let rect = find_empty_rect_containing(&x_pts, j, k, v, grid_size)?;
                let c = compute_u_coeff(&rect, o_ref, oi, oj, i, j, n, grid_size);
                let y = GCGen::from_raw(GCGen::encode(&y_seq), n as u8, 0, 0);
                Some((GCInfGen::new(y, &q * &c), FF2::one()))
            })
        })
        .collect()
}

/// Return a new monomial with the exponents of U_i and U_j exchanged.
fn swap_u(u: &MultiVar<'U', usize>, i: usize, j: usize, n: usize) -> MultiVar<'U', usize> {
    let ei = u.deg_for(i);
    let ej = u.deg_for(j);
    MultiVar::from_iter((0..n).filter_map(|r| {
        let e = if r == i { ej } else if r == j { ei } else { u.deg_for(r) };
        if e > 0 { Some((r, e)) } else { None }
    }))
}

/// Find an empty rectangle from x to y (swap at j, k) that contains v.
/// A rectangle is empty if none of x's points lie in its strict interior.
fn find_empty_rect_containing(
    x_pts: &[GridPoint],
    j: usize,
    k: usize,
    v: GridPoint,
    grid_size: u8,
) -> Option<GridRect> {
    let pj = x_pts[j];
    let pk = x_pts[k];
    for (from, to) in [(pj, pk), (pk, pj)] {
        let rect = GridRect::between(from, to, grid_size);
        if rect.intersects(x_pts, grid_size, true) { continue; }
        if rect.contains(v, grid_size, false) {
            return Some(rect);
        }
    }
    None
}

/// Compute the U-monomial coefficient for a rectangle in the pentagon map.
///
/// The exponent of U_r is:
///   r = i → 0
///   r = j → 1 if Oi is "above" Oj within the rect AND (Oi shifted right) is in the rect
///   other → 1 if the O-mark at column r lies inside the rect
fn compute_u_coeff(
    rect: &GridRect,
    o_pts: &[GridPoint],
    oi: GridPoint,
    oj: GridPoint,
    i: usize,
    j: usize,
    n: usize,
    grid_size: u8,
) -> MultiVar<'U', usize> {
    let l = grid_size as usize;
    let y0 = rect.origin.y as usize;

    MultiVar::from_iter((0..n).filter_map(|r| {
        let e: usize = if r == j {
            let yi = (oi.y as usize + l - y0) % l;
            let yj = (oj.y as usize + l - y0) % l;
            let oi_right = oi.shift(1, 0);
            if yi > yj && rect.contains(oi_right, grid_size, false) { 1 } else { 0 }
        } else if r == i {
            0
        } else {
            if rect.contains(o_pts[r], grid_size, false) { 1 } else { 0 }
        };
        if e > 0 { Some((r, e)) } else { None }
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;
    use crate::{GridComplex, GridDiagram};

    #[test]
    fn test_adjacent_transpositions_identity() {
        let v: Vec<u8> = vec![2, 1, 0];
        assert!(adjacent_transpositions(&v, &v).is_empty());
    }

    #[test]
    fn test_adjacent_transpositions_reverse() {
        let from: Vec<u8> = vec![2, 1, 0];
        let to: Vec<u8> = vec![0, 1, 2];
        let trans = adjacent_transpositions(&from, &to);

        let mut curr = from.clone();
        for (i, j) in &trans {
            curr.swap(*i, *j);
        }
        assert_eq!(curr, to);
    }

    #[test]
    fn test_hom_gen_unknot() {
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCFiltered>::new(gd);
        let z = cx.hom_gen();

        assert!(!z.is_zero());
        assert!(cx.differentiate(&z).is_zero());

        let gh = cx.into_homology(0);
        assert_eq!(gh[0].vectorize(&z).into_vec(), vec![FF2::one()]);
    }

    #[test]
    fn test_hom_gen_trefoil() {
        let gd = GridDiagram::trefoil();
        let cx = GridComplex::<GCFiltered>::new(gd);
        let z = cx.hom_gen();
        
        assert!(!z.is_zero());
        assert!(cx.differentiate(&z).is_zero());

        let gh = cx.into_homology(0);
        assert_eq!(gh[0].vectorize(&z).into_vec(), vec![FF2::one()]);
    }

    #[test]
    fn test_hom_gen_figure8() {
        let gd = GridDiagram::figure8();
        let cx = GridComplex::<GCFiltered>::new(gd);
        let z = cx.hom_gen();
        
        assert!(!z.is_zero());
        assert!(cx.differentiate(&z).is_zero());

        let gh = cx.into_homology(0);
        assert_eq!(gh[0].vectorize(&z).into_vec(), vec![FF2::one()]);
    }
}
