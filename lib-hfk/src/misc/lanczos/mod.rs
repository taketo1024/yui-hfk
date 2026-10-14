//! Block Lanczos over `𝔽₂` (P. Montgomery, *A Block Lanczos Algorithm for Finding
//! Dependencies over GF(2)*, EUROCRYPT 1995), for sparse systems too large for PLUQ.
//!
//! [`find_ng_certificate`] looks for `w` with `aᵀw = 0` and `yᵀw = 1` (so `y ∉ im a`),
//! and [`find_solution`] for `x` with `a·x = y`. Both are randomized, but every vector
//! they return is verified, so a `Some` is always correct and a `None` only means
//! that nothing was found.

mod f2mat;
mod kernel;
mod lanczos;
mod mat64;
mod mix;
mod rng;
#[cfg(test)]
mod test_utils;

pub use lanczos::{find_ng_certificate, find_solution};
