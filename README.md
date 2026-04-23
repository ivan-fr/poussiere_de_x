# Pandrosion Lean 4 Corpus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19712450.svg)](https://doi.org/10.5281/zenodo.19712450)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)
[![Universitas Pandrosion CI](https://github.com/ivan-fr/poussiere_de_x/actions/workflows/ci.yml/badge.svg)](https://github.com/ivan-fr/poussiere_de_x/actions)

Formal verification and numerical illustration of the Pandrosion rational
root-finding map, with a compiled research paper and reproducible figures.

**Cite this work:** DOI [10.5281/zenodo.19712450](https://zenodo.org/records/19712450)

The current primary artifact is:

- [`articles/pandrosion_pth.pdf`](articles/pandrosion_pth.pdf) - distributable compiled paper
- [`latex/pandrosion_pth.tex`](latex/pandrosion_pth.tex) - LaTeX source
- [`lean/Pandrosion.lean`](lean/Pandrosion.lean) - root Lean module importing the corpus

## What This Repository Contains

This project studies the **Pandrosion-Steffensen algorithm** for solving the
complex equation `z^p = x`. The base map is the derivative-free rational
iteration built from the geometric sum `S_p(s) = 1 + s + ... + s^(p-1)`:

```text
h(s) = 1 - (x - 1) / (x * S_p(s))
```

Its fixed points are exactly the `p`-th roots of `1/x`. Composing `h` with the
Aitken-Steffensen `Δ²` accelerator yields a derivative-free iterator of
quadratic order that attains the Kung-Traub efficiency bound for two
evaluations per step.

The repository contains:

- the **`Pandrosion.Core` spine** — 24 Lean 4 modules in `lean/Pandrosion/Core/`
  that build, in dependency order, the algebraic foundations, multi-start
  Voronoï architecture, closed-form linear rate, Kung-Traub optimality,
  complex multiplier, local attraction, almost-everywhere dynamical
  convergence of the Steffensen iterator, the unconditional super-attracting
  Steffensen spine (§§29-34), and the **fully unconditional axiom-clean
  proof** of `RealMcMullenP2` on `ℝ` at `p = 2`, `x > 1` via explicit Möbius
  conjugacy `σ_{2,x} ~ μ²w²`, Böttcher-coordinate squaring-power dichotomy,
  and a formalised Lebesgue-null bad-set lemma (§§35-37);
- the **`Pandrosion.Legacy` companion** — 9 Lean modules in
  `lean/Pandrosion/Legacy/` (non-load-bearing, one-way dependency on `Core`)
  that formalise range invariance of `h` on `[0,1]`, exact Aitken Δ²
  extrapolation on geometric error, the `p = 2` contraction identity and
  the universal `(x−1)/x` contraction bound for all `p ≥ 2`, `HasDerivAt`
  facts for `S_p` and `h_{p,2}`, no-periodic-orbit under a contraction,
  the affine form of the Voronoï bisector, the universal descent
  `D_p < 0`, the cube-root derivatives `P'(r) = −1/5` and `P''(r) = −12/(25r)`,
  the Chebyshev-Halley exclusion, and the anchor-based cube-root variant
  `F_a(s)` with Newton as the degenerate case `a = s`;
- a LaTeX research paper organized around the `Pandrosion.Core` spine plus a
  `Pandrosion.Legacy` companion section;
- Docker Compose tooling for reproducible Lean builds.

## Current Formal Status

As of April 22, 2026:

- Lean toolchain: `leanprover/lean4:v4.7.0`
- Mathlib: `v4.7.0`
- Lean modules under `lean/Pandrosion/`: **34** (24 Core + 9 Legacy + 1 aggregator)
- Total theorems/lemmas: **296** — all audited against the Lean whitelist
- Definitions: **59**
- `lake build Pandrosion` passes; `docker compose run --rm lean-check` passes the max-strict gate with `0 sorry, 0 admit, 0 raw axiom, 0 warning, 0 error, 0 sorryAx, 0 off-whitelist axiom`.
- No executable `sorry` terms; no project-level `axiom` declarations.
- Every declaration's axiom dependency set reduces to `{propext, Classical.choice, Quot.sound}`.

The paper is therefore best read as a formally verified algebraic and
dynamical case study, not as a claim to have resolved any open global problem.
Some figures are numerical or schematic illustrations; the formal claims are
the Lean theorem statements.

## The `Pandrosion.Core` Spine (24 Modules)

The current paper is organized around the 24-module spine in
`lean/Pandrosion/Core/`. Dependency order top-to-bottom:

1. [`Foundations.lean`](lean/Pandrosion/Core/Foundations.lean) — `S_p`, real
   and complex fixed-point theorem `h(s) = s ⇔ s^p = 1/x`, Pandrosion-Steffensen
   accelerator, scaling invariance, multi-start grand master.
2. [`MultiStartBasins.lean`](lean/Pandrosion/Core/MultiStartBasins.lean) —
   Voronoï basins: convexity, closedness, connectedness, bisector frontier;
   constructive McMullen off the boundary.
3. [`QuadraticRate.lean`](lean/Pandrosion/Core/QuadraticRate.lean) — algebraic
   skeleton of Aitken-Steffensen quadratic rate.
4. [`QuadraticComplexity.lean`](lean/Pandrosion/Core/QuadraticComplexity.lean) —
   convergence-to-`ε` bounds from the abstract skeleton.
5. [`VoronoiMeasure.lean`](lean/Pandrosion/Core/VoronoiMeasure.lean) — Voronoï
   boundary has 2D Lebesgue measure zero.
6. [`CyclotomicMcMullen.lean`](lean/Pandrosion/Core/CyclotomicMcMullen.lean) —
   cyclotomic anchor family `γ_k = α · e^(2πik/p)`, unconditional constructive
   McMullen.
7. [`KungTraub.lean`](lean/Pandrosion/Core/KungTraub.lean) — efficiency index,
   Pandrosion attains `KT(2) = √2`, converse optimality (axiom-conditional).
8. [`UniformComplexity.lean`](lean/Pandrosion/Core/UniformComplexity.lean) —
   uniform complexity skeleton across ensembles.
9. [`SuperGrandMaster.lean`](lean/Pandrosion/Core/SuperGrandMaster.lean) —
   abstract super-grand-master for the bounded regime `x > 1`.
10. [`UniformContractionRate.lean`](lean/Pandrosion/Core/UniformContractionRate.lean) —
    closed-form model `λ_p^model(x)`, Taylor residue bound, `p`-uniform complexity.
11. [`ConcreteIteration.lean`](lean/Pandrosion/Core/ConcreteIteration.lean) —
    concrete specialization to Pandrosion-Steffensen sequences.
12. [`MasterAbsolu.lean`](lean/Pandrosion/Core/MasterAbsolu.lean) — top-level
    stitching: five simultaneous conjuncts (anchor distinctness, fixed-point
    property, a.e. Voronoï uniqueness, `p`-uniform complexity, Kung-Traub
    optimality).
13. [`ComplexMultiplier.lean`](lean/Pandrosion/Core/ComplexMultiplier.lean) —
    complex multiplier `h'(α) = 1 - p·α^(p-1)/S_p(α)` in closed form at every
    fixed point.
14. [`LocalAttraction.lean`](lean/Pandrosion/Core/LocalAttraction.lean) — abstract
    local attraction and its quadratic-bound reducer; Steffensen local
    attraction conditional on a local quadratic bound.
15. [`DynamicalConvergence.lean`](lean/Pandrosion/Core/DynamicalConvergence.lean) —
    a.e. dynamical convergence of the Steffensen orbit to a cyclotomic anchor,
    conditional on per-anchor quadratic bounds and a.e. trajectory entry.

**Unconditional Steffensen spine (§§29-34):**

16. [`LinearAsymptotics.lean`](lean/Pandrosion/Core/LinearAsymptotics.lean) —
    `λ_C ≠ 1` at every fixed point; quantitative linear-contraction modulus
    for `h` near `α`.
17. [`QuadraticBound.lean`](lean/Pandrosion/Core/QuadraticBound.lean) —
    closed-form Taylor residue of `h_{p,x}` with explicit `(C_0, r_0)`.
18. [`SteffensenQuadraticBound.lean`](lean/Pandrosion/Core/SteffensenQuadraticBound.lean) —
    **unconditional local quadratic bound on `σ_{p,x}`**; discharges the
    conditional hypothesis from §14/§15.
19. [`SteffensenMcMullenAE.lean`](lean/Pandrosion/Core/SteffensenMcMullenAE.lean) —
    named Prop `McMullenAEEntry`; a.e. solving `z^p = x` modulo McMullen.
20. [`SteffensenExplicitRate.lean`](lean/Pandrosion/Core/SteffensenExplicitRate.lean) —
    computable `(K_α, r_α)`; explicit super-attractive contraction certificate.
21. [`SteffensenGlobalLoglog.lean`](lean/Pandrosion/Core/SteffensenGlobalLoglog.lean) —
    loglog tail inside the basin; a.e. global loglog mod McMullen.

**Real `p = 2` unconditional tower (§§35-37):**

22. [`SteffensenRealP2.lean`](lean/Pandrosion/Core/SteffensenRealP2.lean) —
    Möbius form of `h_{2,x}`; real closed-form multiplier `(1-α)/(1+α) ∈ (0,1)`;
    named Prop `RealMcMullenP2`.
23. [`SteffensenRealMcMullenP2.lean`](lean/Pandrosion/Core/SteffensenRealMcMullenP2.lean) —
    **Möbius conjugacy `σ_{2,x} ~ μ²w²`, cleared form**; explicit rational
    bridge `sigma_p2_explicit`.
24. [`SteffensenRealMcMullenP2Unconditional.lean`](lean/Pandrosion/Core/SteffensenRealMcMullenP2Unconditional.lean) —
    Böttcher coordinate `v`, iterated identity `v(σ^n) = v^(2^n)`, finite
    Julia section on `ℝ`, convergence dichotomy, bad-set countability via
    polynomial fiber-finiteness, Lebesgue-null bad-set lemma, and the
    **fully unconditional axiom-clean `RealMcMullenP2`** theorem
    (`mcmullen_p2_real_unconditional`).

The spine builds with `Compiled: 24 / 24` via `docker compose run --rm
lean-incremental`.

## Repository Layout

```text
.
├── articles/
│   └── pandrosion_pth.pdf            # Current distributable paper
├── latex/
│   ├── pandrosion_pth.tex            # Main paper source
│   ├── pandrosion_pth.pdf            # Local compiled copy
│   └── pandrosion_*.pdf              # Paper figures
├── lean/
│   ├── Pandrosion.lean               # Root import module
│   ├── Pandrosion/
│   │   ├── Core/                     # 24 Lean modules (load-bearing spine)
│   │   ├── Legacy/                   # 9 Lean modules (companion, audited)
│   │   └── Legacy.lean               # Legacy aggregator
│   ├── lakefile.lean
│   ├── lake-manifest.json
│   └── lean-toolchain
├── figures/                          # Additional generated visual assets
├── docker-compose.yml                # Lean build/check services
├── Dockerfile                        # Lean container image
└── tectonic                          # Local Tectonic binary
```

## Reproduce The Lean Build

The easiest path is Docker Compose:

```bash
docker compose run --rm lean-check
```

Expected result for the full corpus (Core + Legacy):

```text
--- Source-tree scan (ALL lean/**/*.lean, comments stripped) ---
  all .lean files (lean/):    37
  Pandrosion/ modules:        34
  .olean built (Pandrosion):  34
  theorems+lemmas:            296
  defs:                       59
  sorry in source:            0
  admit in source:            0
  axiom in source:            0

✅ MAX-STRICT CHECK PASSED — 34 modules, 296 theorems/lemmas, 59 defs, 296 audited
   0 sorry, 0 admit, 0 raw axiom, 0 warning, 0 error, 0 sorryAx, 0 off-whitelist axiom
```

For a clean rebuild with cache cleanup and summary:

```bash
docker compose run --rm lean-build
```

To open an interactive Lean workspace container:

```bash
docker compose run --rm lean
```

The Lean project itself is under `lean/`; inside the container the working
directory is `/workspace/lean`.

## Rebuild The Paper

The repository includes a local `tectonic` binary at the repository root.

```bash
(cd latex && ../tectonic pandrosion_pth.tex)
```

Then copy the compiled PDF to the distributable location:

```bash
cp latex/pandrosion_pth.pdf articles/pandrosion_pth.pdf
```

From inside `latex/`, the equivalent copy command is:

```bash
cp pandrosion_pth.pdf ../articles/pandrosion_pth.pdf
```

## The `Pandrosion.Legacy` Companion (9 Modules)

Alongside the Core spine, `lean/Pandrosion/Legacy/` collects 9 audited
companion modules with strict one-way dependency (`Legacy` imports `Core`,
never the reverse). They add 44 theorems and 7 definitions on top of the
Core corpus and pass the same axiom-audit whitelist. Contents:

- [`Legacy/Basic.lean`](lean/Pandrosion/Legacy/Basic.lean) — range invariance of `h` on `[0,1]` (`h_lt_one`, `h_pos`), orbit-in-`(0,1)` for `p ≥ 2`, exact Aitken Δ² extrapolation (`aitken_perfect_extrapolation`).
- [`Legacy/ContractionP2.lean`](lean/Pandrosion/Legacy/ContractionP2.lean) — the `p = 2` contraction identity `h(s) − h(t) = (x−1)(s−t)/[x(1+s)(1+t)]`, strict distance decrease, monotonicity, uniqueness of the positive fixed point.
- [`Legacy/UniversalContraction.lean`](lean/Pandrosion/Legacy/UniversalContraction.lean) — universal divided difference `Q_p` with `S_p(s) − S_p(t) = (s − t)·Q_p(s, t)`, the `(x−1)/x` contraction bound for every `p ≥ 2` (`contraction_general`), concrete factorisations for `p = 3, 4, 5`.
- [`Legacy/Derivative.lean`](lean/Pandrosion/Legacy/Derivative.lean) — Mathlib-calculus bindings: `HasDerivAt` for `S_p` (all `p`) and for `h_{p,2}`, with closed-form asymptotic rate `h'(s*) = (x−1)/(x(1+s*)²) ∈ (0,1)`.
- [`Legacy/NoCycles.lean`](lean/Pandrosion/Legacy/NoCycles.lean) — no periodic orbit under a contraction (`no_two_cycle`, `no_periodic_orbit`).
- [`Legacy/VoronoiAffine.lean`](lean/Pandrosion/Legacy/VoronoiAffine.lean) — affine coordinate form of the Voronoï bisector; `Fin d` surjection is automatically bijective.
- [`Legacy/Descent.lean`](lean/Pandrosion/Legacy/Descent.lean) — universal descent `D_p = (1/p)·∑ log(cos(kπ/(2p))) < 0` for every `p ≥ 2`, via strict `log cos` bounds on `(0, π/2)`.
- [`Legacy/CubeRoot.lean`](lean/Pandrosion/Legacy/CubeRoot.lean) — cube-root iteration `P(s) = s(s³+4X)/(3s³+2X)`: universal linear rate `P'(r) = −1/5`, quadratic correction `P''(r) = −12/(25r)`, Pandrosion ≠ Newton / Halley / Chebyshev-Halley (polynomial cross-identities, exclusion from the `CH_α` family).
- [`Legacy/AnchorStep.lean`](lean/Pandrosion/Legacy/AnchorStep.lean) — anchor-based cube-root step `F_a(s) = a − (a³−X)/Q(a,s)`, fixed-point at every cube root, Newton as the degenerate case `a = s`, reanchor-at-root, `multistart_step` full-pair fixed-point, `aitken_exact_geometric` variant.

The companion is *not* part of the paper's load-bearing claims but serves as
reusable infrastructure: range bounds, calculus hooks, explicit low-`p`
identities, and the cube-root-specific comparison theorems routinely cited
alongside the Pandrosion iteration.

## Scope Notes

This repository contains several types of evidence:

- Lean theorem proofs: machine-checked formal statements.
- LaTeX exposition: human-readable presentation of those statements and their interpretation.
- Figures: numerical or schematic visualizations, not additional proofs.

**What the paper proves (unconditional):** real and complex fixed-point
characterisation, Steffensen preserves every Pandrosion fixed point, multi-start
grand master, Voronoï basin structure (convex, closed, connected, bisector
frontier), constructive McMullen off the boundary, Voronoï boundary has Lebesgue
measure zero, cyclotomic anchor injectivity, unconditional constructive
McMullen for `z^p = x`, closed-form rate `λ_{p,x} ∈ (0,1)`, complex multiplier
`h'(α) = 1 - p·α^(p-1)/S_p(α)` at every fixed point, Kung-Traub attainment
`E(2,2) = √2`, `p`-uniform complexity bound, Master Absolu, **unconditional
local quadratic bound on the Steffensen iterator** (`SteffensenQuadraticBound`,
§31) and its two corollaries (unconditional local attraction and a.e.
dynamical convergence modulo McMullen), **explicit super-attractive rate**
`K_α · r_α ≤ 1/2` with loglog tail inside the basin, **Möbius conjugacy
`σ_{2,x} ~ μ²w²`** on `ℝ` at `p = 2` in cleared form, **iterated Böttcher
identity** `v(σ^n(s)) = v(s)^(2^n)` and finite Julia section on `ℝ`, the
**convergence dichotomy off the real bad set** (`orbit_enters_basin_off_bad_set`),
**bad-set countability** (`real_bad_set_countable`) via polynomial
fiber-finiteness (`Polynomial.finite_setOf_isRoot` on the quadratic
`A s² + B s + A/x` obtained from the bridge), the corresponding
**Lebesgue-null lemma** (`real_bad_set_measure_zero`), and therefore the
**fully unconditional axiom-clean `RealMcMullenP2(x, α)`** theorem
(`mcmullen_p2_real_unconditional`) for every `x > 1`, `α = 1/√x`, with no
residual hypothesis. All load-bearing theorems audit to the Lean core
whitelist `{propext, Classical.choice, Quot.sound}`.

**What is proven conditional on an explicit axiom:** Kung-Traub converse
optimality (conditional on the Kung-Traub 1974 conjecture, stated as an
explicit `axiom` in `KungTraub.lean`).

**What is proven conditional on explicit hypotheses (no new axioms):**
complex a.e. dynamical convergence at arbitrary `p` — conditional on the
named Prop `McMullenAEEntry p x α` (Fatou-type a.e. trajectory entry into
some local disk). In `ℝ` at `p = 2`, this hypothesis is now **fully and
unconditionally** discharged axiom-clean — including the bad-set
Lebesgue-null premise — via `mcmullen_p2_real_unconditional` (see
`SteffensenRealMcMullenP2Unconditional`, §37).

**What is not claimed:** an unconditional a.e. trajectory-entry theorem in
`ℂ` at arbitrary `p` (the Fatou/McMullen global dynamics beyond the `ℝ`, `p=2`
case); any resolution of classical open problems (Riemann hypothesis, BSD,
Faltings, Smale's 17th, etc.).

The value of the corpus is the exact algebraic identities, the measure-theoretic
coverage theorem, and the dynamical-systems statement — for this particular
rational iteration, in a proof-assistant-checked form.

## Citation

```bibtex
@software{besevic2026pandrosion,
  title  = {Universitas Pandrosion: Formal Verification of a Rational Root-Finding Map and Diophantine Bridges in Lean 4},
  author = {Besevic, Ivan},
  year   = {2026},
  doi    = {10.5281/zenodo.19712450},
  url    = {https://zenodo.org/records/19712450},
  note   = {24-module Pandrosion.Core spine + 9-module Pandrosion.Legacy companion; 296 theorems, 0 sorry, 0 off-whitelist axiom; fully unconditional axiom-clean RealMcMullenP2 on R at p=2, x>1 (bad-set Lebesgue-null lemma formalised via polynomial fiber-finiteness)}
}
```

## License

This work is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

## Author

Ivan Besevic, April 2026 (last updated 2026-04-22).
