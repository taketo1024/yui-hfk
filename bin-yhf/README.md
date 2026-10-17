# yhf

[![crates.io](https://img.shields.io/crates/v/yhf.svg)](https://crates.io/crates/yhf)
[![docs.rs](https://docs.rs/yhf/badge.svg)](https://docs.rs/yhf)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A CLI tool for Heegaard Floer knot homology computations, based on the library
[`yui-hfk`](https://crates.io/crates/yui-hfk).

## Setup

* Install [Rust](https://www.rust-lang.org/tools/install).
* Install `yhf` by

```sh
$ cargo install yhf
```

* Download the grid diagrams of the knots in [KnotInfo](https://knotinfo.org/) by

```sh
$ yhf setup
```

## Usage

```text
Usage: yhf <COMMAND> <KNOT> [OPTIONS]

Commands:
  hfk      the knot Floer homology ĤFK
  gh       the grid homology GH~, GH^ or GH⁻ (-t tilde|hat|minus)
  tau      the Ozsváth–Szabó invariant τ
  upsilon  the Ozsváth–Stipsicz–Szabó invariant Υ
  d        the d-invariants of the p/q-surgery (yhf d <KNOT> <P/Q>)
  g0       the invariant 𝒢₀
  setup    download grid diagrams from KnotInfo into the data directory
```

`<KNOT>` is a knot name such as `3_1`, looked up in the data directory. A grid diagram can
also be given inline as `[os, xs]`, where `os[i]` (resp. `xs[i]`) is the row of the `O`
(resp. `X`) in column `i`. Add `-m` to use the mirror. See [Options](#options) for the other
options, or run `yhf <COMMAND> --help`.

## Examples

### Knot Floer homology `ĤFK`

Command `hfk` computes the ranks of `ĤFK(K)` over `𝔽₂`, bigraded by the Maslov degree `M`
and the Alexander degree `A`. For a grid diagram `G` of `K` with grid number `n`, it is
computed from `GH~(G) ≅ ĤFK(K) ⊗ W^{⊗(n−1)}`, where `W` is 2-dimensional, together with the
symmetry `ĤFK_{i,j}(K) ≅ ĤFK_{i−2j,−j}(K)`.

```sh
$ yhf hfk 3_1
```
```text
 A\M  -2  -1  0 
 1    .   .   1 
 0    .   1   . 
 -1   1   .   . 
```

The same knot can be given by a grid diagram:

```sh
$ yhf hfk "[[2,1,0,4,3],[4,3,2,1,0]]"
```
```text
 A\M  -2  -1  0 
 1    .   .   1 
 0    .   1   . 
 -1   1   .   . 
```

### Grid homology `GH`

Command `gh` computes the grid homology. The default type is `hat`. This command is not
meant for practical use: `GH^` is computed directly, so it takes more time than `ĤFK` by
`hfk`.

```sh
$ yhf gh 3_1
```
```text
GH^ (Hat):
 j\i  -2  -1  0 
 1    .   .   F₂ 
 0    .   F₂  . 
 -1   F₂  .   . 
```

### The `τ`-invariant

Command `tau` computes `τ(K)`.

```sh
$ yhf tau 3_1
1
$ yhf tau 3_1 -m
-1
```

### The `Υ`-invariant

Command `upsilon` computes `Υ_K` from `𝒢₀(K)`, as the piecewise linear function on `[0, 2]`
through the listed vertices `(t, Υ_K(t))`.

```sh
$ yhf upsilon 3_1
PL[(0, 0), (1, -1), (2, 0)]
```

### The `d`-invariants of Dehn surgeries

Command `d` computes the `d`-invariants `d(S³_{p/q}(K), i)` for `0 ≤ i < |p|`, one line
for each `i`. They are given by the invariants `V_k(K)`, which are read from `𝒢₀(K)`. It
computes `τ(K)` first: if `τ(K) ≤ 0` and `R(0, 0)` is realizable, all `V_k(K) = 0`, and
`𝒢₀(K)` is not needed.

```sh
$ yhf d 3_1 3/4
0: -3/2
1: -13/6
2: -13/6
```

### The invariant `𝒢₀`

Command `g0` computes `𝒢₀(K)`. Each line is one region: `R(i, j)` is the quadrant
`(−∞, i] × (−∞, j]`, and `R[(i₁, j₁), (i₂, j₂), …]` is the union of such quadrants.

```sh
$ yhf g0 3_1
R(0, 1)
R(1, 0)
```

## Options

### Known invariants

Known values of invariants make the computation faster.

| option | commands | |
|---|---|---|
| `-g, --g <G>` | `tau`, `upsilon`, `d`, `g0` | The genus. The command aborts if it disagrees with the computed one. `tau` uses it instead of computing the genus. |
| `-G, --g4 <G4>` | `tau`, `upsilon`, `d`, `g0` | The 4-genus. It bounds the candidate regions, and `\|τ\| ≤ g4`. |
| `-t, --tau <TAU>` | `upsilon`, `d`, `g0` | `τ`. It decides the regions `{i ≤ −1} ∪ R(0, m)` for all `m`, and so prunes the candidates. `d` then skips computing `τ`. |

`-G` and `-t` are only checked for consistency with `ĤFK`, so a wrong value can give a wrong
result. `g0` also takes `-R <REGION>` to test only that region, and `--expect <REGION>` to
test the given regions first. The result does not depend on `--expect`.

### Computation

The commands `tau`, `upsilon`, `d` and `g0` test candidate regions by solving linear systems
over `𝔽₂`. These options control the solver and show its progress.

- **`--log <N>`** — `0` (off, the default), `1` (info), `2` (debug) or `3` (trace), written to
  stderr. At `1`, it shows the input diagram, as `[os, xs]` and as a grid, with its grid
  number and writhe. Then it shows the steps: `ĤFK` and the genus, the candidate regions,
  the numbers of generators by bidegree, and the result of each test: `OK` (realizable),
  `NG` (not realizable) or `RESIGN` (not decided). At `2`, it also shows the sizes of the
  linear systems and the progress of the solver.

  ```text
  $ yhf g0 3_1 --log 1
  [INFO  yui_hfk::grid_diagram::diagram] diagram: [[2, 1, 0, 4, 3], [4, 3, 2, 1, 0]]
       j\i  0  1  2  3  4
       4    X  .  .  O  .
       3    .  X  .  .  O
       2    O  .  X  .  .
       1    .  O  .  X  .
       0    .  .  O  .  X

  [INFO  yui_hfk::grid_diagram::diagram]   grid number: 5
  [INFO  yui_hfk::grid_diagram::diagram]   writhe: 3
  [INFO  yui_hfk::g_zero::calc] compute HFK ...
  ...
  [INFO  yui_hfk::g_zero::calc] try: R(0, 1) (remain: 0)
  [INFO  yui_hfk::g_zero::calc] OK: R(0, 1) (the last candidate)
  ```

  The invariants of a knot and its mirror (`-m`) determine each other, for example
  `τ(K*) = −τ(K)`, and `𝒢₀(K*)` is the dual of `𝒢₀(K)`. But the cost on the two sides can
  differ a lot, and the numbers of generators help to choose the side.
- **`--max-size <N>`** — a region whose linear system has dimension over `N` (default
  `1000000`) is not tested, but resigned. If the result needs the region, the command fails
  with `Resigned`. `tau` first tries the mirror. A larger `N` decides more regions, at the
  cost of time and memory.
- **`--try-fast-ng`** — before the full system, tries the smaller system obtained by
  specializing the `U`-variables. It can show that a region is not realizable (`NG`), also
  when the full system is over `--max-size` (if the smaller one is not). It never shows that
  a region is realizable.
- **`--mode <SOLVER>`** — the solver for the linear systems.
  `lanczos` (the default) is block Lanczos: it uses little memory, and is randomized but
  checks its answers. `lanczos(n)` tries up to `n` random starts.
  `pluq` is a deterministic elimination, which may fill in on large systems. It solves them
  in steps of `100000`, or of `n` with `pluq(n)`.

## Data directory

Grid diagrams are stored as `<data_dir>/grid_diagrams/<name>.json`, where `<data_dir>` is

1. `$YUI_DATA_DIR` if set.
2. Otherwise the platform user-data dir:
    - macOS: `~/Library/Application Support/yui/`
    - Linux: `${XDG_DATA_HOME:-~/.local/share}/yui/`
    - Windows: `%APPDATA%\yui\`

`yhf setup` keeps existing files unless `--force` is given. `--out DIR` writes to another
directory, and `--in FILE` reads a local copy of the KnotInfo database (`.xls` or
`.xls.zip`) instead of downloading it.

## License

`yhf` is released under the [MIT license](https://opensource.org/licenses/MIT).

---

*This README was generated by [Claude](https://www.anthropic.com/claude).*
