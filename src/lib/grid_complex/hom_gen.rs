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

        let mut z: Lc<GCInfGen, FF2> = Lc::from(GCInfGen::new(MultiVar::<'U', usize>::one(), x0));

        curr_os = o0;
        for (i, j) in trans {
            z = apply_pentagon(z, i, j, &curr_os, grid_size);
            curr_os.swap(i, j);
        }

        // Replace dummy generators with actual generators from the diagram
        z.into_iter().filter_map(|(inf_gen, coeff)| {
            self.gens.gen_from_code(inf_gen.x.code).map(|actual| {
                (GCInfGen::new(inf_gen.u, actual.clone()), coeff)
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
            let x_seq = inf_gen.x.sequence();
            let x_pts = inf_gen.x.points();
            let q = swap_u(&inf_gen.u, i, j, n);
            let o_ref = o_pts.as_slice();
            (0..n).filter_map(move |k| {
                if k == j { return None; }
                let mut y_seq = x_seq.clone();
                y_seq.swap(j, k);
                let rect = find_empty_rect_containing(&x_pts, j, k, v, grid_size)?;
                let c = compute_u_coeff(&rect, o_ref, oi, oj, i, j, n, grid_size);
                let y = GCGen::from_raw(GCGen::encode(&y_seq), n as u8, 0, 0);
                Some((GCInfGen::new(&q * &c, y), FF2::one()))
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
    use itertools::Itertools;

    use super::*;
    use num_traits::Zero;
    use crate::{GridComplex, GridDiagram};

    // All grid diagrams of grid number n, skipping the degenerate ones that put
    // an O and an X in the same cell.
    fn all_diagrams(n: u8) -> impl Iterator<Item = GridDiagram> {
        let perms = (0..n).permutations(n as usize).collect_vec();
        perms.clone().into_iter().cartesian_product(perms)
            .filter(|(os, xs)| os.iter().zip(xs).all(|(o, x)| o != x))
            .map(|(os, xs)| GridDiagram::from_sequences(&os, &xs))
    }

    // The defining property of hom_gen: it is a cycle. Run exhaustively over
    // every grid diagram up to size 5 — this is what exercises the pentagons
    // that wrap around the torus, which the named knots below never produce.
    #[test]
    fn test_hom_gen_is_cycle_exhaustive() {
        for n in 3..=5u8 {
            let mut count = 0;
            for gd in all_diagrams(n) {
                let cx = GridComplex::<GCFiltered>::new(gd);
                let z = cx.hom_gen();

                assert!(!z.is_zero(), "n = {n}: hom_gen vanished");
                assert!(cx.differentiate(&z).is_zero(), "n = {n}: d(z) != 0");
                count += 1;
            }
            assert!(count > 0, "n = {n}: no diagrams generated");
        }
    }

    #[test]
    fn test_hom_gen_is_cycle_for_named_knots() {
        for name in ["0_1", "3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3"] {
            let gd = GridDiagram::test_data(name).unwrap();
            let cx = GridComplex::<GCFiltered>::new(gd);
            let z = cx.hom_gen();

            assert!(!z.is_zero(), "{name}: hom_gen vanished");
            assert!(cx.differentiate(&z).is_zero(), "{name}: d(z) != 0");
        }
    }

    #[test]
    fn test_hom_gen_is_cycle_for_mirrors() {
        for name in ["3_1", "4_1", "5_1", "5_2"] {
            let gd = GridDiagram::test_data(name).unwrap().rotate90();
            let cx = GridComplex::<GCFiltered>::new(gd);
            let z = cx.hom_gen();

            assert!(!z.is_zero(), "{name} mirror: hom_gen vanished");
            assert!(cx.differentiate(&z).is_zero(), "{name} mirror: d(z) != 0");
        }
    }

    // hom_gen really does emit monomials that are not powers of U₀ alone, so the
    // distinction between alg_deg (= -a₀) and the total degree is observable
    // here. This is the witness that made the algebraic-degree bug visible.
    #[test]
    fn test_hom_gen_emits_mixed_monomials() {
        let cx = GridComplex::<GCFiltered>::new(GridDiagram::figure8());
        let z = cx.hom_gen();

        let mixed = z.iter().filter(|(g, _)| {
            g.u.deg_for(0) == 0 && g.u.total_deg() > 0
        }).count();

        assert!(mixed > 0, "expected a term with no U₀ factor but some Uᵢ factor");

        // Such a term has alg_deg = 0 while -total_deg < 0.
        assert!(z.iter().any(|(g, _)| g.alg_deg() != -(g.u.total_deg() as isize)));
    }

    #[test]
    fn test_swap_u() {
        let u = MultiVar::<'U', usize>::from([1usize, 2, 3]);
        let v = swap_u(&u, 0, 2, 3);
        assert_eq!(v.deg_for(0), 3);
        assert_eq!(v.deg_for(1), 2);
        assert_eq!(v.deg_for(2), 1);
        // swapping twice is the identity
        assert_eq!(swap_u(&v, 0, 2, 3), u);
    }

    #[test]
    fn test_swap_u_fixes_untouched_vars() {
        let u = MultiVar::<'U', usize>::from([0usize, 5, 0, 7]);
        let v = swap_u(&u, 0, 2, 4);
        // U₀ and U₂ are both absent, so the monomial is unchanged
        assert_eq!(v, u);
    }

    #[test]
    fn test_adjacent_transpositions_sorts_random_perms() {
        for n in 2..=5u8 {
            for from in (0..n).permutations(n as usize) {
                for to in (0..n).permutations(n as usize) {
                    let trans = adjacent_transpositions(&from, &to);
                    let mut curr = from.clone();
                    for (i, j) in &trans {
                        assert_eq!(*j, i + 1, "must be adjacent");
                        curr.swap(*i, *j);
                    }
                    assert_eq!(curr, to);
                }
            }
        }
    }

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
        assert_eq!(gh[0].vectorize(&z).into_dense(), vec![FF2::one()]);
    }

    #[test]
    fn test_hom_gen_trefoil() {
        let gd = GridDiagram::trefoil();
        let cx = GridComplex::<GCFiltered>::new(gd);
        let z = cx.hom_gen();
        
        assert!(!z.is_zero());
        assert!(cx.differentiate(&z).is_zero());

        let gh = cx.into_homology(0);
        assert_eq!(gh[0].vectorize(&z).into_dense(), vec![FF2::one()]);
    }

    #[test]
    fn test_hom_gen_figure8() {
        let gd = GridDiagram::figure8();
        let cx = GridComplex::<GCFiltered>::new(gd);
        let z = cx.hom_gen();
        
        assert!(!z.is_zero());
        assert!(cx.differentiate(&z).is_zero());

        let gh = cx.into_homology(0);
        assert_eq!(gh[0].vectorize(&z).into_dense(), vec![FF2::one()]);
    }
}
