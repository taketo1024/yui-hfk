# Review: `hom_gen.rs` vs SS2019 Section 5.3

Date: 2026-08-06
Scope: `src/lib/grid_complex/hom_gen.rs` checked against `references/SS2019.pdf`
(Sections 4.1, 4.2, 5.2, 5.3, 5.4), plus the geometric primitives
(`grid_diagram/rect.rs`, `point.rs`, `diagram.rs`) and the downstream consumers
(`g_zero/calc.rs`, `complex_ext.rs`, `region.rs`, `tau.rs`) that `hom_gen` feeds.

## Verdict

**The algorithm in `hom_gen.rs` is correct.** The pentagon-count combinatorics
were re-derived from scratch on the torus (degenerating `γ → β`, checking the
one-quadrant corner condition at both crossing points of `β ∩ γ`); the code's
rule matches the derivation exactly.

**A separate bug was found in the surrounding `g_zero` machinery** — the
algebraic degree was computed from the total `U`-degree instead of the
`U₀`-exponent alone. It is outside `hom_gen.rs` but sits under everything
`hom_gen` feeds. See the section below; it has been fixed.

## Verified, step by step

### 1. The starting cycle `x₀` (Prop 5.4)

`hom_gen.rs:21-28` places the O's at `o0[i] = n−1−i` and takes `x₀` at the
top-left corners of the O's. The paper's anti-diagonal is `row = n−i−2`, so the
code's configuration is the paper's **translated up by one row** on the torus.
Harmless — the complex is translation-invariant. Prop 5.4's proof was checked
directly for the code's configuration: for every transposition the two
connecting rectangles are both empty and O-free, so they cancel over `𝔽₂` and
`∂x₀ = 0`. The hardcoded `(0, 0)` gradings on the dummy generators are never
consulted; the final `filter_map` swaps in the real generators. ✓

### 2. The transposition sequence

`adjacent_transpositions` is a correct bubble sort; each intermediate
O-configuration is a permutation, so each step is a valid column commutation.
Sorting `o0 → o1` directly is equivalent to the paper's "transform G to G₀,
then reverse". Ignoring the X's is justified exactly as the paper says: neither
`∂` nor `Φ` sees them. ✓

### 3. Pentagon = rectangle containing `v`

With `(i, j) = (c, c+1)`, the circle `β` is `x = 2j`, and the code models the
`β ∩ γ` crossing as `v = (2j, oj.y)`. Derivation: with `γ` based left of `β`
and bulging around the two O's, both crossing points lie in the interior of
row `b` (the right O's row, pre-move). Since rectangle edges sit at even
heights, any point in row b's interior is combinatorially equivalent — so one
`v` at height `2b+1` correctly represents the crossing. A pentagon from x to y
exists iff the connecting rectangle has row b in its row range, on either side
of `β`. Three facts make the code exact here:

- `contains(·, ·, false)` is **closed on all edges** (`rect.rs:69-75`), so `v`
  (whose x-coordinate lies on a vertical edge) is admitted on both the
  left-of-β and right-of-β rectangles;
- the two rectangles of `Rect(x, y)` have **disjoint open row ranges**, so `v`
  lands in exactly one — no missed and no double-counted pentagon;
- emptiness is checked on strict interiors, which agrees with
  pentagon-emptiness because the ε-slivers contain no even-coordinate points. ✓

### 4. The U-exponents (`compute_u_coeff`)

Deriving the pentagon-vs-rectangle O-count difference: the bulge sliver around
the right O (target label `U_i`) is always excised — so `ε_i = 0`, matching
`r == i → 0`. The left O (target label `U_j`) is picked up by the γ-arc sliver
iff its height lies in the cyclic interval from `2b+1` up to the rectangle's
top edge. The code's condition `yi > yj && rect.contains(oi_right, …)` — both
heights measured cyclically from `rect.origin.y` — is **exactly** this interval
condition, uniformly for both left and right rectangles. All other O's: plain
containment, and markings (odd coords) never touch edges, so the closed-edge
convention is immaterial for them. ✓

**Remark.** The paper's phrasing ("O_i is located at the left-upper of
O_{i+1}…") is planar and becomes ambiguous when Π wraps around the torus. The
code's cyclic measurement from the rectangle's bottom edge is the correct
toroidal refinement — the code is more precise than the paper's prose here.
Figure 5's caption could note this.

### 5. `swap_u`

`Φ` is module-equivariant under the relabeling that follows the physical
markings, which swaps `U_i ↔ U_j`. Applying `swap_u` to the incoming monomial
and multiplying by the target-labeled pentagon coefficient is right. ✓

### 6. Discarded coefficients

`apply_pentagon` drops the incoming `FF2` coefficient. Valid: yui-1.0's `Lc` is
always reduced, so stored coefficients are nonzero = 1. ✓

### 7. Empirical corroboration

`∂z = 0` and `[z] = generator` pass for unknot / trefoil / figure-8, and the τ
tests (`3_1`, `4_1`, `5_1`, `5_2` + mirrors) exercise the full chain with knots
where τ takes distinct values.

## Bug found and fixed: the algebraic degree (outside `hom_gen.rs`)

Raised here first as an open question; confirmed by the author as a **porting
error** from the original Swift implementation, and fixed in commit `be6b793`
on branch `fix-algebraic-degree`.

### The error

The paper (Sections 4.2 and 5.2) defines

```
Alg (U₁^{a₁} ⋯ Uₙ^{aₙ} x) = −a₁           ← the U₁-exponent only
Alex(U₁^{a₁} ⋯ Uₙ^{aₙ} x) = A(x) − Σᵢ aᵢ  ← the total degree
```

`Alg` uses only `a₁` because `C⁻(G)` is an `𝔽[U]`-module via `U = U₁`. The
Swift original agrees: `algebraicDegree = -left.exponent[0]`. The Rust port
had `alg_deg = −Σᵢ aᵢ`, the total degree.

Origin: commit `91a860e` ("add GCIngGen::alg_deg") introduced `alg_deg` and, in
the same edit, rewrote `alex_deg` from `A(x) − Σaᵢ` to `A(x) + alg_deg()`.
`alex_deg` kept returning the right number, but only because it was expressed
through an `alg_deg` that happened to equal it — two distinct quantities
collapsed into one.

Independent confirmation that this was a bug and not a convention: `shift_z`
in `g_zero/calc.rs` multiplies by `U₀^s`, i.e. the code already treats `U₀` as
the distinguished variable, which is inconsistent with `alg = −Σaᵢ`.

### Blast radius

`alg_deg` is not merely an accessor. The whole `(alg, alex)` indexing in
`g_zero/complex_ext.rs` assumed `alg = −total`, which makes the total degree
**determined** by the point (`e = −i`). Under the correct definition it is not:
a point `(alg, alex) = (i, j)` at Maslov degree `k` draws from **every**
`e ≥ −i`, with raw generators at `(k + 2e, j + e)` and multiplicity

```
#{monomials in U₁ … U_{n−1} of degree e + i} = C(n − 2 + e + i, e + i).
```

Rewritten accordingly: `complementary_points`, `count_inf_gens_in`,
`inf_gens_in`, `count_raw_gens_in`, `raw_gens_in`, plus a new `a₀`-fixed
generator `inf_gens_e_at` in `grid_complex/complex.rs`. The enumeration is now
driven from the raw bidegrees rather than the raw Alexander range, since
correct inflated points can fall below that range.

### Direction and observed effect

Since `−Σaᵢ ≤ −a₀` and regions are down-closed, `C_R^{buggy} ⊇ C_R^{correct}`:
the old region test was too permissive and `Q_R` too small, so the code could
**accept regions that are not realizable** — `𝒢₀` too small and `τ` too small,
never too large. Concretely, with corner `(−1, ∞)` the old convention made
`alg = 0` mean "no `U` at all", so `Q_R` was just the raw generators; correct
is `a₀ = 0` with `U₁ … U_{n−1}` free, which is far larger. This matches the
paper's remark that `A` "becomes extremely large as `n` increases".

The divergence is real, not theoretical: for the figure-8, `hom_gen` produces a
term with `a₀ = 0` but `U₁`-degree > 0, which the two conventions place at
different algebraic degrees.

**However, computed results are unchanged for every knot in the test data.**
`τ` and `𝒢₀` were computed before and after the fix for
`3_1, 4_1, 5_1, 5_2, 6_1, 6_2, 6_3` and agree exactly
(`τ = 1, 0, 2, 1, 0, 1, 0`, all matching known values). The bug needs mixed
monomials at the critical degrees, which grids this small do not produce.
Exhibiting a difference requires larger knots — especially ones with rich
`CFK^∞` such as `8_19` — so a full recomputation is still warranted.

### Regression guards

Five tests now discriminate the two conventions; all five fail if `alg_deg` is
reverted to the total degree:
`chain::test_alg_degree`, `chain::test_alg_degree_no_unit`,
`chain::test_support_grid`, `gc_gen::test_inf_gen_degrees`,
`gc_gen::test_alg_deg_ignores_higher_vars`.

## Minor notes

- `hom_gen`'s final `filter_map` silently drops terms absent from `self.gens`.
  For the current callers this is **provably sound**: `build_complex`'s
  predicate evaluates each generator at exactly the inflation degree that
  `z`-terms and the `b`-vector use, so dropped terms lie inside every candidate
  region and project to zero in every `Q_R`, and `contains_chain` answers are
  unaffected. But the contract is invisible at the call site; a doc comment on
  `hom_gen` ("caller must guarantee dropped generators are irrelevant, see
  `build_complex`") would prevent future misuse.
- `z.is_zero()` in `from_gd` is only logged. It happens to be sound
  (all-terms-dropped implies every candidate contains z and is genuinely
  realizable), but a comment saying so would stop it looking like a swallowed
  error.
- Suggested test: a property test asserting `∂(hom_gen()) = 0` over random
  grid diagrams (random O/X permutations, n ≤ 7). The torus-wrapping pentagons
  are the only part not exercised by the three curated knots, and that is
  exactly where a regression would hide.
