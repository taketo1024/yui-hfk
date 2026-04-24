//! YUI-HFK: Heegaard Floer Knot Homology in Rust
//!
//! A Rust implementation of the grid-based Heegaard Floer knot homology,
//! ported from the Swift library `swm-hfk`.
//!
//! See arXiv:2002.09210 for mathematical background.

pub mod hfk;
pub mod grid_diagram;
pub mod grid_complex;
pub mod g_zero;
pub mod misc;

pub use hfk::HFK;
pub use grid_diagram::{GridDiagram, GridPoint, GridRect};
pub use grid_complex::{GCGen, GridComplex, GCType, GCTilde, GCHat, GCMinus, GCFiltered};
