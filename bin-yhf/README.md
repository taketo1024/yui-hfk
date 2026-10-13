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
  hfk    the knot Floer homology ĤFK
  gh     the grid homology GH~, GH^ or GH⁻ (-t tilde|hat|minus)
  tau    the Ozsváth–Szabó invariant τ
  g0     the invariant 𝒢₀
  setup  download grid diagrams from KnotInfo into the data directory
```

`<KNOT>` is a knot name such as `3_1`, looked up in the data directory. A grid diagram can
also be given inline as `[os, xs]`, where `os[i]` (resp. `xs[i]`) is the row of the `O`
(resp. `X`) in column `i`. Add `-m` to use the mirror, and `--log 1` (up to `3`) to follow
the computation. Run `yhf <COMMAND> --help` for all options.

## Examples

Command `hfk` computes the ranks of `ĤFK(K)` over `𝔽₂`, bigraded by the Maslov degree `M`
and the Alexander degree `A`.

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

Command `gh` computes the grid homology. The default type is `hat`.

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

Command `tau` computes `τ(K)`.

```sh
$ yhf tau 3_1
1
$ yhf tau 3_1 -m
-1
```

Command `g0` computes `𝒢₀(K)`. Each line is one region: `R(i, j)` is the quadrant
`(−∞, i] × (−∞, j]`, and `R[(i₁, j₁), (i₂, j₂), …]` is the union of such quadrants.

```sh
$ yhf g0 3_1
R(0, 1)
R(1, 0)
```

For `tau` and `g0`, the option `-G/--g4` supplies a known value of the 4-genus, and `g0`
also accepts a known value of `τ` by `-t/--tau`. These prune the candidate regions. `--max-size` caps the size of the linear systems that are attempted.

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
