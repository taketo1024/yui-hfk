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

Use [`yhf`](https://crates.io/crates/yhf) to use the program from the command line.

## What it computes

- `ĤFK(K)` — the knot Floer homology, bigraded by Maslov and Alexander degree.
- `g(K)` — the Seifert genus, read off from the top Alexander degree.
- `τ(K)` — the Ozsváth–Szabó concordance invariant.
- `Υ_K` — the Ozsváth–Stipsicz–Szabó invariant \[OSS17\], a piecewise linear function on `[0, 2]`.
- `d(S³_{p/q}(K))` — the `d`-invariants of the Dehn surgeries along `K` \[OS03a\].
- `𝒢₀(K)` — the concordance invariant of \[Sat23\], a finite set of closed regions in `ℤ²`
  from which `τ(K)`, `Υ_K` and the `d`-invariants of Dehn surgeries can be read.

## Key types

### Grid diagrams

- **`GridDiagram`** — an `n × n` grid with one `O` and one `X` in each row and column.
  Built by `from_sequences(os, xs)`, where `os[i]` (resp. `xs[i]`) is the row of the `O`
  (resp. `X`) in column `i`, or by `from_arc_presentation`. `load(name)` reads a diagram
  by knot name (see [Grid diagram data](#grid-diagram-data)), and `rotate90()` gives a
  diagram of the mirror. `sequences()` gives `[os, xs]` back, `to_table_string()` draws the
  grid, and `log_diagram()` logs both, with the grid number and the writhe.
- **`GridPoint`**, **`GridRect`** — points and rectangles on the grid torus.

### Grid complexes

- **`GridComplex<V>`** — the grid complex of a `GridDiagram` over `𝔽₂`, where the type
  parameter `V` selects the variant:
  `GCTilde` (`GC~`, \[Book\] Def 4.4.1), `GCHat` (`GC^`, Def 4.6.12),
  `GCMinus` (`GC⁻`, Def 4.6.1) and `GCFiltered` (the filtered complex, Def 13.2.1).
  `into_complex()` and `into_homology()` give the chain complex and its homology, and
  `log_generators()` logs the numbers of generators by bidegree.
- **`GCGen`** — a grid state, with its Maslov and Alexander degrees.
  **`GCInfGen`** — a grid state times a monomial in the `U`-variables.

### Invariants

- **`HFK`** — the bigraded ranks of `ĤFK(K)` over `𝔽₂`, built by `HFK::from(&gd)`.
  It is computed from `GH~(G) ≅ ĤFK(K) ⊗ W^{⊗(n−1)}`, where `W` is 2-dimensional,
  together with the symmetry `ĤFK_{i,j}(K) ≅ ĤFK_{i−2j,−j}(K)`.
- **`genus(&gd)`** — the Seifert genus `g(K)`, the top Alexander degree of `GH~(G)`.
- **`tau(&gd)`** — the Ozsváth–Szabó invariant `τ(K)`, by the regions `{i ≤ −1} ∪ R(0, m)`
  (\[SS20\] Thm 2.20). If `τ ≤ −1` or the search resigns, it is computed on the mirror, by
  `τ(K*) = −τ(K)`. `tau_with_config` also accepts the genus and the 4-genus when they are
  known, and a `G0CalcConfig`.
- **`upsilon(&gd)`** — `Υ_K` as an `Upsilon`, a piecewise linear function with rational
  vertices, read from `𝒢₀(K)` (\[SS20\] Thm 2.20). `Upsilon::from_g0` reads it from a known `𝒢₀(K)`.
- **`d(&gd, p, q)`** — the `d`-invariants `d(S³_{p/q}(K), i)` for `0 ≤ i < |p|`, from the
  `V_k` of \[NW15\]. `v_seq(&gd)` gives `[V_0, V_1, …]`: it computes `τ` first, and all
  `V_k = 0` if `τ ≤ 0` and `R(0, 0)` is realizable. Otherwise they are read from `𝒢₀(K)`
  (\[SS20\] Thm 2.20). `v_seq_from_g0` reads them from a known `𝒢₀(K)`, and
  `v_seq_from_regions` tests the regions `R(m, k+m)` directly, for validation.
- **`G0Calc`** — computes `𝒢₀(K)` as a list of `Region`s. A `Region` is a finite union of
  quadrants `(−∞, i] × (−∞, j]`, stored as an antichain of `Corner`s. `G0CalcConfig`
  bounds the size of the linear systems. Regions that could not be decided within
  these bounds are reported by `resigned()`, and `tau`, `upsilon` and `d` then return
  `Err(G0CalcError::Resigned)`.

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
`yhf hfk 3_1`, `yhf tau 3_1`, `yhf upsilon 3_1`, `yhf d 3_1 5/2`, `yhf g0 3_1`, etc.
See its [README](https://github.com/taketo1024/yui-hfk/blob/main/bin-yhf/README.md) for the
commands and examples.

## Tests

```bash
cargo test                                                          # unit tests and small knots
cargo test --release --test thick_11 -- --ignored --nocapture       # slow: one file at a time
cargo test --release --test thick_12_t23 tau_all -- --ignored --nocapture
```

The ignored tests check `τ` (`tau_*`) and `𝒢₀` (`g0_*`), with the genus, on knots that are
homologically thick and non-slice, or have `τ = 0`:

- [`tests/thick_11.rs`](https://github.com/taketo1024/yui-hfk/blob/main/lib-hfk/tests/thick_11.rs)
  and [`tests/tau0_11.rs`](https://github.com/taketo1024/yui-hfk/blob/main/lib-hfk/tests/tau0_11.rs):
  the 44 knots of Table 1 in \[SS20\] (grid number 7–11), against the original Swift
  implementation.
- `tests/thick_12_t*.rs` (one file for each `𝒢₀`-type), `tests/tau0_12.rs` and
  `tests/thick_12_heavy.rs`: the 199 prime knots with 12 crossings, against the `𝒢₀`-types
  of \[SS20\] (v2), with `g`, `g4` and `|τ|` from KnotInfo.

The `τ` tests are cheaper, about 1h45m for all of them on an Apple M1 Ultra. The `𝒢₀` tests
take much longer, up to about a day for one knot of `thick_12_heavy.rs` (see the `ignore`
reasons there), so run them one file, or one test, at a time. Expected values and grid
diagrams are inlined, so the tests need no external data.

## References

- \[Book\] P. Ozsváth, A. Stipsicz, Z. Szabó,
  *Grid Homology for Knots and Links*, AMS Math. Surveys and Monographs 208 (2015).
- \[MOST07\] C. Manolescu, P. Ozsváth, Z. Szabó, D. Thurston,
  *On combinatorial link Floer homology*, Geom. Topol. 11 (2007).
- \[NW15\] Y. Ni, Z. Wu, *Cosmetic surgeries on knots in S³*, J. Reine Angew. Math. 706 (2015).
- \[OS03a\] P. Ozsváth, Z. Szabó, *Absolutely graded Floer homologies and intersection forms
  for four-manifolds with boundary*, Adv. Math. 173 (2003).
- \[OSS17\] P. Ozsváth, A. Stipsicz, Z. Szabó, *Concordance homomorphisms from knot Floer
  homology*, Adv. Math. 315 (2017).
- \[Sat23\] K. Sato, *The ν⁺-equivalence classes of genus one knots*, Selecta Math. (N.S.) 29
  (2023).
- \[SS20\] T. Sano, K. Sato, *Computation of the Υ-invariant and the d-invariants of Dehn
  surgeries for knots with up to 12 crossings*, [arXiv:2002.09210](https://arxiv.org/abs/2002.09210v2).

## License

`yui-hfk` is released under the [MIT license](https://opensource.org/licenses/MIT).

---

*This README was generated by [Claude](https://www.anthropic.com/claude).*
