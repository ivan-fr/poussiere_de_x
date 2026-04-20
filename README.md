# Pandrosion Lean 4 Corpus

Formal verification and numerical illustration of the Pandrosion rational
root-finding map, with a compiled research paper and reproducible figures.

The current primary artifact is:

- [`articles/pandrosion_paper.pdf`](articles/pandrosion_paper.pdf) - distributable compiled paper
- [`latex/pandrosion_paper.tex`](latex/pandrosion_paper.tex) - LaTeX source
- [`lean/Pandrosion.lean`](lean/Pandrosion.lean) - root Lean module importing the corpus

## What This Repository Contains

This project studies the Pandrosion iteration for root extraction. For the
general `p`-th root setting the core map is written using the geometric sum
`S_p(s)`:

```text
h(s) = 1 - (x - 1) / (x * S_p(s))
```

For the cubic specialization, the main rational map is:

```text
P_X(s) = s * (s^3 + 4X) / (3s^3 + 2X)
```

The repository contains:

- a Lean 4 formalization of algebraic, dynamical, Diophantine, matrix, and spectral identities;
- a LaTeX research paper describing the formal corpus;
- generated figures, including a final ten-figure proof gallery;
- Python scripts for numerical and visual exploration;
- Docker Compose tooling for reproducible Lean builds.

## Current Formal Status

As of April 20, 2026:

- Lean toolchain: `leanprover/lean4:v4.7.0`
- Mathlib: `v4.7.0`
- Lean modules under `lean/Pandrosion/`: `100`
- Top-level `theorem` declarations under `lean/Pandrosion/` plus root import file: `678`
- `lake build Pandrosion` passes through Docker Compose.
- No executable `sorry` terms were found in the corpus; the only `sorry` occurrence is in prose inside a comment.
- Two explicit axioms remain in [`lean/Pandrosion/SpectralLimit.lean`](lean/Pandrosion/SpectralLimit.lean), used for spectral-limit statements:
  - `integral_log_cos_eq`
  - `D_eq_closed`

The paper is therefore best read as a formally verified algebraic and
dynamical case study, not as a claim to have resolved any open global problem.
Some figures are numerical or schematic illustrations; the formal claims are
the Lean theorem statements.

## Proof Highlights

The final section of the paper contains a ten-figure proof gallery generated
from [`scripts/generate_proof_gallery_figures.py`](scripts/generate_proof_gallery_figures.py):

1. Universal local derivative `P'(r) = -1/5`
2. Chebyshev-Halley family exclusion
3. Kinematic residual conservation
4. Exclusion of periodic orbits under contraction
5. Pell-Pandrosion integer norm amplification
6. Cross-determinant separation of consecutive approximants
7. Voronoi basin stability under contraction
8. Hermitian preservation and spectral confinement
9. DFT character orthogonality
10. Effective irrationality bound and norm explosion

The generated gallery files live in `latex/fig_proof_gallery_*.pdf` and are
included by `latex/pandrosion_paper.tex`.

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
│   ├── Pandrosion/                   # 100 Lean modules
│   ├── lakefile.lean
│   ├── lake-manifest.json
│   └── lean-toolchain
├── scripts/
│   └── generate_proof_gallery_figures.py
├── verification/                     # Python verification/exploration scripts
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

Expected result:

```text
Compiled: 100 / 100
ALL 100 MODULES OK
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

## Regenerate The Proof Gallery Figures

Requirements:

- Python 3
- `numpy`
- `matplotlib`

On this machine, `/usr/bin/python3` has the needed scientific packages. The
Matplotlib cache is pointed at `/tmp` to avoid user-cache permission issues:

```bash
MPLCONFIGDIR=/tmp/mplconfig /usr/bin/python3 scripts/generate_proof_gallery_figures.py
```

Generic command if your default Python environment has the dependencies:

```bash
python3 scripts/generate_proof_gallery_figures.py
```

After regenerating figures, rebuild the paper:

```bash
(cd latex && ../tectonic pandrosion_paper.tex)
```

## Useful Lean Modules

Some central modules:

- [`lean/Pandrosion/HalleyComparison.lean`](lean/Pandrosion/HalleyComparison.lean) - Newton/Halley comparison and universal derivative `-1/5`
- [`lean/Pandrosion/ChebyshevHalleyExclusion.lean`](lean/Pandrosion/ChebyshevHalleyExclusion.lean) - exclusion from the Chebyshev-Halley family
- [`lean/Pandrosion/ResidualConservation.lean`](lean/Pandrosion/ResidualConservation.lean) - kinematic residual conservation
- [`lean/Pandrosion/NoCycles.lean`](lean/Pandrosion/NoCycles.lean) - no finite cycles under contraction
- [`lean/Pandrosion/ThueBridge.lean`](lean/Pandrosion/ThueBridge.lean) - norm amplification and cross-determinants
- [`lean/Pandrosion/VoronoiInvariance.lean`](lean/Pandrosion/VoronoiInvariance.lean) - Voronoi convexity and basin stability
- [`lean/Pandrosion/HermitianPreservation.lean`](lean/Pandrosion/HermitianPreservation.lean) - Hermitian preservation for matrix products
- [`lean/Pandrosion/DFTDecomposition.lean`](lean/Pandrosion/DFTDecomposition.lean) - roots-of-unity cancellation and DFT identities
- [`lean/Pandrosion/EffectiveIrrationality.lean`](lean/Pandrosion/EffectiveIrrationality.lean) - effective Liouville-type lower bound

## Python Verification And Exploration

The `verification/` directory contains exploratory numerical scripts and
stress tests used while developing the paper. These are not the formal proof
source; Lean is the authoritative formal layer.

Examples:

```bash
python3 verification/verification_complex.py
python3 verification/verification_optimality.py
python3 verification/final_stress_test.py
```

Some scripts may require scientific Python packages such as `numpy`,
`matplotlib`, or `scipy`.

## Scope Notes

This repository contains several types of evidence:

- Lean theorem proofs: machine-checked formal statements.
- LaTeX exposition: human-readable presentation of those statements and their interpretation.
- Figures: numerical or schematic visualizations, not additional proofs.
- Python verification scripts: exploratory tests and numerical experiments.

The project does not claim to solve Smale's 17th problem, the Riemann
hypothesis, abc, Roth, or other global open problems. The value of the corpus
is in the exact algebraic identities and certified structural bridges that it
formalizes for this particular rational iteration.

## Citation

```bibtex
@article{besevic2026pandrosion,
  title  = {The Pandrosion-Steffensen Iteration: Formal Verification of a Rational Root-Finding Map in Lean 4},
  author = {Besevic, Ivan},
  year   = {2026},
  note   = {Lean 4 formalization and preprint}
}
```

## License

This work is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

## Author

Ivan Besevic, April 2026
