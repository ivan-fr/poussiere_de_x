# Universitas Pandrosion

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19733297.svg)](https://doi.org/10.5281/zenodo.19733297)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)
[![Universitas Pandrosion CI](https://github.com/ivan-fr/poussiere_de_x/actions/workflows/ci.yml/badge.svg)](https://github.com/ivan-fr/poussiere_de_x/actions)

Formal verification (Lean 4) of the Pandrosion-Steffensen iteration for
`z^p = x`, with a companion series of nine research papers developing the
adaptive-anchor generalisation to arbitrary complex polynomials and to
multivariate polynomial systems (the proper setting of Smale's 17th problem),
plus reproducible numerical experiments.

**Cite this work:** DOI [10.5281/zenodo.19733297](https://zenodo.org/records/19733297)

## The Nine-Paper Corpus (`articles/`)

| # | File | Scope |
|---|------|-------|
| 0 | [`0pandrosion_pth.pdf`](articles/0pandrosion_pth.pdf) | **Pth**: Pandrosion-Steffensen for `z^p = x`, formally verified in Lean 4 (the primary load-bearing artifact). |
| 1 | [`1pandrosion_smale.pdf`](articles/1pandrosion_smale.pdf) | **Part I**: Generalised Pandrosion operator `P_{P,z_0}(z)` for any univariate polynomial `P ∈ C[z]`; multivariate extension `F: C^n → C^n` via the Schmidt slope matrix (Smale 17 setting, with full attribution to Lairez 2017 and Beltrán–Pardo 2009). |
| 2 | [`2pandrosion_smale.pdf`](articles/2pandrosion_smale.pdf) | **Part II**: Smale's Mean Value Conjecture reformulated as a Pandrosion ratio bound, with exact proofs for `d = 2, 3`. |
| 3 | [`3pandrosion_smale.pdf`](articles/3pandrosion_smale.pdf) | **Part III**: The (classical Lagrange–Sylvester) vanishing identity `Σ 1/P'(α_k) = 0` and its implications for Smale MVC. |
| 4 | [`4pandrosion_smale.pdf`](articles/4pandrosion_smale.pdf) | **Part IV**: Pandrosion inverse, fiber vanishing identity, proof of Smale MVC for `d = 3` (already classical, Smale 1981). |
| 5 | [`5pandrosion_smale.pdf`](articles/5pandrosion_smale.pdf) | **Part V**: Smale MVC for `d = 4` via a resultant inequality — the conjecture for `d = 4` was previously established (Tischler 1989, Beardon–Minda–Ng 2002); we provide an alternative path via Pandrosion reduction. |
| 6 | [`6pandrosion_smale.pdf`](articles/6pandrosion_smale.pdf) | **Part VI**: Smale MVC for `d = 5` via a compactness argument (strict local maximum at the extremal `z^5 + z` + Vieta domination at infinity + Lojasiewicz on the discriminant locus). Supported by a sympy symbolic Hessian certificate and 70 001 numerical tests with 0 violations. **Subject to peer review**: if the compactness argument holds, this would be a new result for `d = 5`. |
| 7 | [`7pandrosion_smale.pdf`](articles/7pandrosion_smale.pdf) | **Part VII**: Derivative-free adaptive Pandrosion scheme for univariate root-finding — breaking McMullen's barrier via a non-holomorphic fallback (Armijo backtracking and analytic γ-damped variants). |
| 8 | [`8pandrosion_smale.pdf`](articles/8pandrosion_smale.pdf) | **Part VIII**: The geometric-decreasing start strategy eliminates the empirical worst-case phase transition for Pandrosion-`T_2`. |

### Paper 0 in detail (Lean-verified core)

This project studies the **Pandrosion-Steffensen algorithm** for solving the
complex equation `z^p = x`. The base map is the derivative-free rational
iteration built from the geometric sum `S_p(s) = 1 + s + ... + s^(p-1)`:

```text
h(s) = 1 - (x - 1) / (x * S_p(s))
```

Its fixed points are exactly the `p`-th roots of `1/x`. Composing `h` with the
Aitken-Steffensen `Δ²` accelerator yields a derivative-free iterator of
quadratic order that attains the Kung-Traub efficiency bound for two
evaluations per step. The formal verification of this paper in Lean 4 is the
primary load-bearing artifact of the repository.

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
- a LaTeX research paper (Paper 0) organised around the `Pandrosion.Core` spine
  plus a `Pandrosion.Legacy` companion section;
- a companion series of eight LaTeX papers (Parts I–VIII) extending the scheme
  to arbitrary complex polynomials and polynomial systems;
- reproducible Python scripts for every numerical experiment cited in the
  papers (regression fits, Hessian certificates, multivariate benchmarks,
  Smale MVC sweeps);
- Docker Compose tooling for reproducible Lean builds.

## Current Formal Status

As of 2026-04-24:

- Lean toolchain: `leanprover/lean4:v4.7.0`
- Mathlib: `v4.7.0`
- Lean modules under `lean/Pandrosion/`: **34** (24 Core + 9 Legacy + 1 aggregator)
- Total theorems/lemmas: **296** — all audited against the Lean whitelist
- Definitions: **59**
- `lake build Pandrosion` passes; `docker compose run --rm lean-check` passes the max-strict gate with `0 sorry, 0 admit, 0 raw axiom, 0 warning, 0 error, 0 sorryAx, 0 off-whitelist axiom`.
- No executable `sorry` terms; no project-level `axiom` declarations.
- Every declaration's axiom dependency set reduces to `{propext, Classical.choice, Quot.sound}`.
- LaTeX corpus: **9 papers** (Paper 0 formally verified in Lean; Parts I–VIII are companion analytical and numerical work, not formally verified).

The Lean corpus is best read as a formally verified algebraic and dynamical
case study on `z^p = x`. The Parts I–VIII companion papers are
*not* formally verified and should be read as traditional mathematical exposition
with explicit numerical certification.

## The `Pandrosion.Core` Spine (24 Modules)

Paper 0 is organised around the 24-module spine in
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

## Companion Papers (Parts I–VIII): Scope and Non-Claims

The eight companion papers develop the adaptive-anchor Pandrosion framework
beyond monomial `z^p = x` to arbitrary univariate polynomials (Part I) and
to polynomial systems `F : C^n → C^n` (Part I §3.5, the actual setting of
Smale's 17th problem). They are analytical and numerical work, not formally
verified in Lean.

### Smale's 17th Problem: what we do and do not claim

Smale's 17th problem (1998) asks for a polynomial-time algorithm to compute
approximate zeros of polynomial systems `F : C^n → C^n` in average. It was
**resolved unconditionally** by:

- Beltrán and Pardo (2009, *J. Amer. Math. Soc.* **22**): probabilistic polynomial
  average time via `μ`-theory.
- Bürgisser and Cucker (2011, *Annals of Mathematics* **174**): deterministic
  `N^{O(log log N)}`.
- Lairez (2017, *Found. Comput. Math.* **17**): fully deterministic average
  polynomial time.

**Our contribution** is a derivative-free algorithmic alternative based on
the multivariate Pandrosion (Schmidt slope matrix) operator, with:

- explicit per-orbit constants and a unified framework that specialises to
  Newton (dynamic anchor) and to Steffensen-type schemes (fixed anchor);
- empirical convergence on Kostlan–Smale random systems up to `(n, d) = (5, 4)`
  (Bézout degree `D = 1024`) with 100% success rate on tested cells;
- a non-holomorphic Armijo fallback mechanism (Part VII) that circumvents
  McMullen's topological barrier for purely rational iterators;
- a geometric-decreasing multistart strategy (Part VIII) that empirically
  eliminates the worst-case phase transition at `d = 128` for i.i.d. Gaussian
  inputs.

We do **not** claim to improve upon the complexity bounds of Lairez or
Beltrán–Pardo; the corresponding basin-entry bound at the `D → ∞` scale
remains an open question in our framework. The algorithmic contribution
stands independently.

### Smale's Mean Value Conjecture (MVC, 1981): what we do

Smale's MVC (1981) is a *different* conjecture from Problem 17 — it is a
univariate inequality on critical points of polynomials. The cases:

- `d = 2, 3`: classical, Smale (1981).
- `d = 4`: proved by Tischler (1989) and others.
- `d = 5`: **open** in the published literature. Part VI of this corpus
  presents a compactness-based proof candidate with: (i) a symbolic Hessian
  certificate (sympy), (ii) 70 001 numerical tests with 0 violations, and
  (iii) a Lojasiewicz argument for the discriminant locus. This result is
  presented as work-in-progress subject to peer review by a specialist (e.g.
  David Minda, Aimo Hinkkanen, or Gerald Schmieder).

Parts II–V provide alternative Pandrosion-language formulations of the
already-classical cases `d ≤ 4`, with explicit credit to the original
references.

## Repository Layout

```text
.
├── articles/                           # Nine distributable compiled papers
│   ├── 0pandrosion_pth.pdf             # Paper 0 (Lean-verified, Pth)
│   ├── 1pandrosion_smale.pdf           # Part I: univariate + multivariate Pandrosion
│   ├── 2pandrosion_smale.pdf           # Part II: Smale MVC reformulation
│   ├── 3pandrosion_smale.pdf           # Part III: Lagrange-Sylvester vanishing identity
│   ├── 4pandrosion_smale.pdf           # Part IV: Pandrosion inverse, fiber identity
│   ├── 5pandrosion_smale.pdf           # Part V: Smale MVC for d=4
│   ├── 6pandrosion_smale.pdf           # Part VI: Smale MVC for d=5 (candidate)
│   ├── 7pandrosion_smale.pdf           # Part VII: Armijo + γ-damped fallback
│   └── 8pandrosion_smale.pdf           # Part VIII: geometric-decreasing starts
├── latex/                              # LaTeX sources + Python experiments
│   ├── 0pandrosion_pth.tex             # Paper 0 source (alias of pandrosion_pth.tex)
│   ├── pandrosion_pth.tex              # Paper 0 source (primary name, Lean-verified)
│   ├── 1pandrosion_smale.tex           # Part I source
│   ├── ...                             # Parts II–VIII sources
│   ├── benchmark_smale17_v4.py         # Multivariate Pandrosion benchmark (n=1..4, d=2..5)
│   ├── benchmark_multivariate_extended.py # Extended benchmark up to (n,d)=(5,4), D=1024
│   ├── explore_higher_order.py         # T_n hierarchy (Richardson) comparison
│   ├── explore_start_strategies.py     # Strategy A (Cauchy) vs B (GeomDecr) vs others
│   ├── explore_mvc_d5.py               # MVC d=5 numerical certificate (70 001 tests)
│   ├── prove_d5_hessian.py             # Symbolic Hessian for MVC d=5 (sympy)
│   ├── prove_armijo_O1.py              # Armijo amortised j_accept distribution
│   ├── measure_worst_case_freq.py      # Worst-case % on UNI random polynomials
│   ├── explore_gamma_damped*.py        # Analytic γ-damped descent step
│   └── *.pdf                           # Compiled copies + figures
├── lean/
│   ├── Pandrosion.lean                 # Root import module
│   ├── Pandrosion/
│   │   ├── Core/                       # 24 Lean modules (load-bearing spine)
│   │   ├── Legacy/                     # 9 Lean modules (companion, audited)
│   │   └── Legacy.lean                 # Legacy aggregator
│   ├── lakefile.lean
│   ├── lake-manifest.json
│   └── lean-toolchain
├── docker-compose.yml                  # Lean build/check services
├── Dockerfile                          # Lean container image
└── tectonic                            # Local Tectonic binary
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

## Rebuild The Nine Papers

The repository includes a local `tectonic` binary at the repository root.

To rebuild all nine papers from `latex/` and copy them to `articles/`:

```bash
cd latex
for f in 0pandrosion_pth 1pandrosion_smale 2pandrosion_smale 3pandrosion_smale \
         4pandrosion_smale 5pandrosion_smale 6pandrosion_smale 7pandrosion_smale \
         8pandrosion_smale; do
  ../tectonic "$f.tex"
  cp "$f.pdf" "../articles/$f.pdf"
done
```

To rebuild a single paper (e.g. Paper 0):

```bash
(cd latex && ../tectonic pandrosion_pth.tex)
cp latex/pandrosion_pth.pdf articles/0pandrosion_pth.pdf
```

## Reproduce The Numerical Experiments

Each `explore_*.py`, `benchmark_*.py`, `prove_*.py`, and
`measure_*.py` script in `latex/` is self-contained. They depend on
`numpy` and `scipy` (and `sympy` for `prove_d5_hessian.py`). Create a
local environment and run:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy sympy
cd latex
python3 benchmark_multivariate_extended.py  # (n,d) up to (5,4), D=1024
python3 explore_mvc_d5.py                   # Smale MVC d=5: 70 001 tests
python3 prove_d5_hessian.py                 # Symbolic Hessian certificate
python3 prove_armijo_O1.py                  # Armijo amortised cost
# ... etc.
```

Each script prints a self-contained summary and is referenced by exact name
from the corresponding paper's text.

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

The companion is *not* part of Paper 0's load-bearing claims but serves as
reusable infrastructure: range bounds, calculus hooks, explicit low-`p`
identities, and the cube-root-specific comparison theorems routinely cited
alongside the Pandrosion iteration.

## Scope Notes

This repository contains several types of evidence:

- **Lean theorem proofs (Paper 0)**: machine-checked formal statements.
- **LaTeX exposition (Paper 0 + Parts I–VIII)**: human-readable presentation of the statements and their interpretation. Only Paper 0 is formally verified.
- **Python numerical experiments**: reproducible computations referenced by paper name and script name.
- **Figures**: numerical or schematic visualisations, not additional proofs.

**What Paper 0 proves (unconditional):** real and complex fixed-point
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

**What is not claimed in Paper 0:** an unconditional a.e. trajectory-entry
theorem in `ℂ` at arbitrary `p` (the Fatou/McMullen global dynamics beyond
the `ℝ`, `p=2` case); any resolution of classical open problems (Riemann
hypothesis, BSD, Faltings, etc.).

**What Parts I–VIII contribute (not formally verified):** a derivative-free
algorithmic framework for polynomial root-finding with an empirical complexity
bound, explicitly complementary (not competitive) to Lairez's deterministic
average-case resolution of Smale 17. Parts II–V reformulate the classical
cases of Smale MVC (`d ≤ 4`) in Pandrosion language. Part VI presents a
proof candidate for Smale MVC at `d = 5` (the first case open in the
literature), supported by a symbolic Hessian certificate and exhaustive
numerical verification, subject to peer review.

The value of the corpus is the exact algebraic identities, the measure-theoretic
coverage theorem, the dynamical-systems statements (Paper 0, Lean-verified),
the algorithmic framework with empirical certification (Parts I, VII, VIII),
and the Smale MVC exposition with a new proof candidate at `d = 5` (Part VI).

## Citation

```bibtex
@software{besevic2026pandrosion,
  title  = {Universitas Pandrosion: Formally Verified Pandrosion-Steffensen Iteration (Lean 4) and Companion Papers on Polynomial Root-Finding and Smale's Mean Value Conjecture},
  author = {Besevic, Ivan},
  year   = {2026},
  doi    = {10.5281/zenodo.19733297},
  url    = {https://zenodo.org/records/19733297},
  note   = {Nine-paper corpus. Paper 0: 34-module Lean 4 verification of Pandrosion-Steffensen for z^p = x, 296 theorems, 0 sorry, 0 off-whitelist axiom, fully unconditional axiom-clean RealMcMullenP2 on R at p=2, x>1. Parts I-VIII: analytical and numerical companions on the adaptive-anchor generalisation, multivariate Smale 17 setting (complementing Lairez 2017), non-holomorphic fallback, and Smale MVC (including a d=5 proof candidate in Part VI subject to peer review).}
}
```

## License

This work is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

## Author

Ivan Besevic, April 2026 (last updated 2026-04-24).
