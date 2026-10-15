# yui-hfk

[![crates.io](https://img.shields.io/crates/v/yui-hfk.svg)](https://crates.io/crates/yui-hfk)
[![docs.rs](https://docs.rs/yui-hfk/badge.svg)](https://docs.rs/yui-hfk)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`yui-hfk` computes *Heegaard Floer knot homology* from grid diagrams, written in
[Rust](https://www.rust-lang.org). It is built on the [`yui`](https://github.com/taketo1024/yui)
libraries (`yui-core`, `yui-homology`, `yui-matrix`).

Grid homology is a combinatorial description of knot Floer homology, introduced by
Manolescu, Ozsváth, Sarkar, Szabó and Thurston. From a grid diagram `G` of a knot `K`
one builds a chain complex whose homology recovers `HFK(K)`, with no analysis involved.

## What it computes

- `ĤFK(K)` — the knot Floer homology, bigraded by Maslov and Alexander degree.
- `GH~`, `GH^`, `GH⁻` — the grid homology variants.
- `g(K)` — the Seifert genus, read off from the top Alexander degree.
- `τ(K)` — the Ozsváth–Szabó concordance invariant.
- `𝒢₀(K)` — the concordance invariant of \[Sat23\], a finite set of closed regions in `ℤ²`
  from which `τ(K)`, `Υ_K` and the `d`-invariants of Dehn surgeries can be read.

## Key types

### Grid diagrams

- **`GridDiagram`** — an `n × n` grid with one `O` and one `X` in each row and column.
  Built by `from_sequences(os, xs)`, where `os[i]` (resp. `xs[i]`) is the row of the `O`
  (resp. `X`) in column `i`, or by `from_arc_presentation`. `load(name)` reads a diagram
  by knot name (see [Grid diagram data](#grid-diagram-data)), and `rotate90()` gives a
  diagram of the mirror.
- **`GridPoint`**, **`GridRect`** — points and rectangles on the grid torus.

### Grid complexes

- **`GridComplex<V>`** — the grid complex of a `GridDiagram` over `𝔽₂`, where the type
  parameter `V` selects the variant:
  `GCTilde` (`GC~`, \[Book\] Def 4.4.1), `GCHat` (`GC^`, Def 4.6.12),
  `GCMinus` (`GC⁻`, Def 4.6.1) and `GCFiltered` (the filtered complex, Def 13.2.1).
  `into_complex()` and `into_homology()` give the chain complex and its homology.
- **`GCGen`** — a grid state, with its Maslov and Alexander degrees.
  **`GCInfGen`** — a grid state times a monomial in the `U`-variables.

### Invariants

- **`HFK`** — the bigraded ranks of `ĤFK(K)` over `𝔽₂`, built by `HFK::from(&gd)`.
  It is computed from `GH~(G) ≅ ĤFK(K) ⊗ W^{⊗(n−1)}`, where `W` is 2-dimensional,
  together with the symmetry `ĤFK_{i,j}(K) ≅ ĤFK_{i−2j,−j}(K)`.
- **`genus(&gd)`** — the Seifert genus `g(K)`, the top Alexander degree of `GH~(G)`.
- **`tau(&gd)`** — the Ozsváth–Szabó invariant `τ(K)`. `tau_with_config` also accepts the
  genus and the 4-genus when they are known, and a `G0CalcConfig`.
- **`G0Calc`** — computes `𝒢₀(K)` as a list of `Region`s. A `Region` is a finite union of
  quadrants `(−∞, i] × (−∞, j]`, stored as an antichain of `Corner`s. `G0CalcConfig`
  bounds the size of the linear systems. Regions that could not be decided within
  these bounds are reported by `resigned()`, and `tau` then returns `Err(G0CalcError::Resigned)`.

## Example

```rust
use yui_hfk::{GridDiagram, HFK, tau};

let gd = GridDiagram::from_sequences(&[2, 1, 0, 4, 3], &[4, 3, 2, 1, 0]);
assert_eq!(tau(&gd), Ok(1));
assert_eq!(HFK::from(&gd).genus(), 1);
```

## Grid diagram data

`GridDiagram::load(name)` reads from an external user-data directory:

1. `$YUI_DATA_DIR` if set.
2. Otherwise the platform user-data dir:
    - macOS: `~/Library/Application Support/yui/`
    - Linux: `${XDG_DATA_HOME:-~/.local/share}/yui/`
    - Windows: `%APPDATA%\yui\`

Grid diagrams live in `<data_dir>/grid_diagrams/<name>.json`, each holding `[os, xs]`
as in `GridDiagram::from_sequences`. To populate the directory from the
[KnotInfo](https://knotinfo.org/) database, use the CLI [`yhf`](https://crates.io/crates/yhf):

```sh
cargo install yhf
yhf setup
```

## CLI

A command-line front-end ships in the sibling crate [`yhf`](https://crates.io/crates/yhf):
`yhf hfk 3_1`, `yhf tau 3_1`, `yhf g0 3_1`, etc.

## Tests

```bash
cargo test                                       # unit tests
cargo test --release -- --ignored --nocapture    # slow: G_0 and tau over knots up to 11 crossings
```

The ignored tests in [`tests/thick_11.rs`](https://github.com/taketo1024/yui-hfk/blob/main/lib-hfk/tests/thick_11.rs) and [`tests/tau0_11.rs`](https://github.com/taketo1024/yui-hfk/blob/main/lib-hfk/tests/tau0_11.rs) check `𝒢₀`, `τ` and the genus
against results from the original Swift implementation, for the 44 knots of Table 1 in \[SS\]
(grid number 7–11).
Expected values and grid diagrams are inlined, so they need no external data.

## References

- \[MOST07\] C. Manolescu, P. Ozsváth, Z. Szabó, D. Thurston,
  *On combinatorial link Floer homology*, Geom. Topol. 11 (2007).
- \[Book\] P. Ozsváth, A. Stipsicz, Z. Szabó,
  *Grid Homology for Knots and Links*, AMS Math. Surveys and Monographs 208 (2015).
- \[Sat23\] K. Sato, *A full-twist inequality for the ν⁺-invariant*, and the `𝒢₀` invariant.
- \[SS\] T. Sano, K. Sato, *An algorithm for computing the Υ-invariant and the
  d-invariants of Dehn surgeries*, [arXiv:2002.09210](https://arxiv.org/abs/2002.09210).

## License

`yui-hfk` is released under the [MIT license](https://opensource.org/licenses/MIT).

---

*This README was generated by [Claude](https://www.anthropic.com/claude).*
