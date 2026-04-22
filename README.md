# Pandrosion Lean 4 Corpus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19689482.svg)](https://doi.org/10.5281/zenodo.19689482)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)
[![Universitas Pandrosion CI](https://github.com/ivan-fr/poussiere_de_x/actions/workflows/ci.yml/badge.svg)](https://github.com/ivan-fr/poussiere_de_x/actions)

Formal verification and numerical illustration of the Pandrosion rational
root-finding map, with a compiled research paper and reproducible figures.

**Cite this work:** DOI [10.5281/zenodo.19689482](https://zenodo.org/records/19689482)

The current primary artifact is:

- [`articles/pandrosion_paper.pdf`](articles/pandrosion_paper.pdf) - distributable compiled paper
- [`latex/pandrosion_paper.tex`](latex/pandrosion_paper.tex) - LaTeX source
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
  Steffensen spine (§§29-34), and the axiom-clean discharge of
  `RealMcMullenP2` on `ℝ` at `p = 2`, `x > 1` via explicit Möbius conjugacy
  `σ_{2,x} ~ μ²w²` and Böttcher-coordinate squaring-power dichotomy (§§35-37);
- a broader **legacy corpus** of further Lean modules under `lean/Pandrosion/`
  exploring algebraic, Diophantine, matrix, and spectral identities;
- a LaTeX research paper organized around the `Pandrosion.Core` spine;
- Docker Compose tooling for reproducible Lean builds.

## Current Formal Status

As of April 20, 2026:

- Lean toolchain: `leanprover/lean4:v4.7.0`
- Mathlib: `v4.7.0`
- Lean modules under `lean/Pandrosion/`: `125`
- Top-level `theorem` declarations under `lean/Pandrosion/` plus root import file: `740+`
- `lake build Pandrosion` passes through Docker Compose.
- No executable `sorry` terms were found in the corpus; the only `sorry` occurrence is in prose inside a comment.
- No project-level `axiom` declarations remain under `lean/Pandrosion/`.
- [`lean/Pandrosion/SpectralLimit.lean`](lean/Pandrosion/SpectralLimit.lean) keeps two classical analytic inputs as explicit Lean variables, used for spectral-limit statements:
  - `integral_log_cos_eq`
  - `D_eq_closed`
  This makes the dependency honest and local: the convergence proof is formalized relative to the closed-form identity, while the identity itself is treated as a standard analytic input rather than a project-declared axiom.

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
    Julia section on `ℝ`, convergence dichotomy, **unconditional axiom-clean
    `RealMcMullenP2` modulo bad-set Lebesgue-null lemma**.

The spine builds with `Compiled: 24 / 24` via `docker compose run --rm
lean-incremental`.

## Repository Layout

```text
.
├── articles/
│   └── pandrosion_paper.pdf          # Current distributable paper
├── latex/
│   ├── pandrosion_paper.tex          # Main paper source
│   ├── pandrosion_paper.pdf          # Local compiled copy
│   ├── pandrosion_master.tex/.pdf    # Larger working manuscript
│   └── fig_*.pdf                     # Paper figures
├── lean/
│   ├── Pandrosion.lean               # Root import module
│   ├── Pandrosion/                   # 125 Lean modules
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

Expected result for the `Pandrosion.Core` spine:

```text
Compiled: 24 / 24
✅ INCREMENTAL OK — 24 modules compiled
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
(cd latex && ../tectonic pandrosion_paper.tex)
```

Then copy the compiled PDF to the distributable location:

```bash
cp latex/pandrosion_paper.pdf articles/pandrosion_paper.pdf
```

From inside `latex/`, the equivalent copy command is:

```bash
cp pandrosion_paper.pdf ../articles/pandrosion_paper.pdf
```

## Useful Lean Modules (Legacy Corpus)

Beyond the `Pandrosion.Core` spine, the legacy corpus under `lean/Pandrosion/`
contains further exploratory modules on the cubic specialization `P_X(s) =
s(s^3+4X)/(3s^3+2X)` and on Diophantine / matrix / spectral identities. These
are *not* included in the current paper's scope but remain in the repository:

- [`lean/Pandrosion/HalleyComparison.lean`](lean/Pandrosion/HalleyComparison.lean) - Newton/Halley comparison and universal derivative `-1/5`
- [`lean/Pandrosion/ChebyshevHalleyExclusion.lean`](lean/Pandrosion/ChebyshevHalleyExclusion.lean) - exclusion from the Chebyshev-Halley family
- [`lean/Pandrosion/ResidualConservation.lean`](lean/Pandrosion/ResidualConservation.lean) - kinematic residual conservation
- [`lean/Pandrosion/NoCycles.lean`](lean/Pandrosion/NoCycles.lean) - no finite cycles under contraction
- [`lean/Pandrosion/ThueBridge.lean`](lean/Pandrosion/ThueBridge.lean) - norm amplification and cross-determinants
- [`lean/Pandrosion/VoronoiInvariance.lean`](lean/Pandrosion/VoronoiInvariance.lean) - Voronoi convexity and basin stability
- [`lean/Pandrosion/HermitianPreservation.lean`](lean/Pandrosion/HermitianPreservation.lean) - Hermitian preservation for matrix products
- [`lean/Pandrosion/DFTDecomposition.lean`](lean/Pandrosion/DFTDecomposition.lean) - roots-of-unity cancellation and DFT identities
- [`lean/Pandrosion/EffectiveIrrationality.lean`](lean/Pandrosion/EffectiveIrrationality.lean) - effective Liouville-type lower bound

## Legacy Frontier Modules (April 2026)

The following modules are part of the broader legacy corpus, not the current
paper's scope. They remain in the repository as exploratory material:

- [`lean/Pandrosion/LittlewoodSimultaneous.lean`](lean/Pandrosion/LittlewoodSimultaneous.lean) - Littlewood frontier `liminf n·‖nα‖·‖nβ‖ = 0` via Schmidt subspace + dim-2 Vojta
- [`lean/Pandrosion/BealOrbital.lean`](lean/Pandrosion/BealOrbital.lean) - Beal / Fermat-Catalan frontier (coprime `A^x + B^y = C^z` with exponents ≥ 3) via abc global frontier + Catalan orbital + Bilu-Tichý
- [`lean/Pandrosion/OppenheimErgodic.lean`](lean/Pandrosion/OppenheimErgodic.lean) - Oppenheim density of indefinite quadratic form values via Riemann gyroscopic attractor + Szemerédi ergodic resonance + Voronoi global density
- [`lean/Pandrosion/LeopoldtPadic.lean`](lean/Pandrosion/LeopoldtPadic.lean) - Leopoldt non-vanishing of the p-adic regulator via p-adic Hensel + non-Archimedean orbital + Baker linear forms
- [`lean/Pandrosion/VojtaMainFrontier.lean`](lean/Pandrosion/VojtaMainFrontier.lean) - Vojta Main Conjecture frontier unifying Roth, abc, Schmidt, Faltings, and Lang under one orbital roof
- [`lean/Pandrosion/ArtinPrimitiveRoot.lean`](lean/Pandrosion/ArtinPrimitiveRoot.lean) - Artin primitive-root infinitude (GRH-conditional) via Riemann gyroscopic + DFT character orthogonality
- [`lean/Pandrosion/LehmerTotient.lean`](lean/Pandrosion/LehmerTotient.lean) - Lehmer totient frontier (`φ(n) | n-1 ⇒ n` prime) via Lehmer spectral limit + Smyth orbital
- [`lean/Pandrosion/GoldfeldAverageRank.lean`](lean/Pandrosion/GoldfeldAverageRank.lean) - Goldfeld average-rank 1/2 for quadratic twists via BSD attractor rank + Brauer-Siegel
- [`lean/Pandrosion/TateAlgebraicCycles.lean`](lean/Pandrosion/TateAlgebraicCycles.lean) - Tate cycle-rank / Galois-invariant rank equality via effective Faltings + matrix Diophantine
- [`lean/Pandrosion/MasonStothersFLT.lean`](lean/Pandrosion/MasonStothersFLT.lean) - Mason-Stothers polynomial abc and function-field FLT via abc global frontier + Thue bridge

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
identity** `v(σ^n(s)) = v(s)^(2^n)` and finite Julia section on `ℝ`, and the
**convergence dichotomy off the real bad set** (`orbit_enters_basin_off_bad_set`)
that discharges `RealMcMullenP2(x, α)` for `x > 1`, `α = 1/√x` modulo a
standalone Lebesgue-null lemma for a countably-constructed bad set. All
load-bearing theorems audit to the Lean core whitelist `{propext,
Classical.choice, Quot.sound}`.

**What is proven conditional on an explicit axiom:** Kung-Traub converse
optimality (conditional on the Kung-Traub 1974 conjecture, stated as an
explicit `axiom` in `KungTraub.lean`).

**What is proven conditional on explicit hypotheses (no new axioms):**
complex a.e. dynamical convergence at arbitrary `p` — conditional on the
named Prop `McMullenAEEntry p x α` (Fatou-type a.e. trajectory entry into
some local disk). In `ℝ` at `p = 2`, this hypothesis is discharged axiom-clean
(see `SteffensenRealMcMullenP2Unconditional`, §37).

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
  doi    = {10.5281/zenodo.19689482},
  url    = {https://zenodo.org/records/19689482},
  note   = {24-module Pandrosion.Core spine, zero sorry; unconditional RealMcMullenP2 discharge on R at p=2}
}
```

## License

This work is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

## Author

Ivan Besevic, April 2026
