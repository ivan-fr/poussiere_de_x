# Universitas Pandrosion

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19757311.svg)](https://doi.org/10.5281/zenodo.19757311)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Universitas Pandrosion CI](https://github.com/ivan-fr/poussiere_de_x/actions/workflows/ci.yml/badge.svg)](https://github.com/ivan-fr/poussiere_de_x/actions)

Universitas Pandrosion is a research repository for finite-slope Pandrosion
methods: derivative-free root extraction, multivariate polynomial systems,
dynamic atlases, Riemann-chart normalization, and exploratory extensions toward
physics and zeta geometry.

The repository contains three kinds of work:

- Lean 4 formalization under [`lean/`](lean/);
- standalone NumPy engines under [`flow/`](flow/);
- LaTeX research papers, PDFs, and generated figures under [`latex/`](latex/)
  and [`scripts/`](scripts/).

Only the Lean artifacts are formal proofs.  The Python engines and LaTeX papers
are reproducible research code and mathematical exposition unless a result is
explicitly proved in the relevant source.

## Current Contents

### Research papers

The current paper series is stored in [`latex/tex/`](latex/tex/) with compiled
PDFs in [`latex/pdf/`](latex/pdf/).  It includes papers `000` through `028`.

Main lines:

- `000`: monomial Pandrosion, contraction, Thales scaling, Riemann inversion,
  and Steffensen acceleration;
- `001` to `019`: finite-slope Halley, anchored Pandrosion, dynamic atlases,
  tensor/inverse-jet geometry, geodesic probe flows, and discrete convergence;
- `020` to `023`: Smale 17 style Pandrosion atlas programs;
- `024`: analog Pandrosion root extraction with Thales-Riemann normalization;
- `025` and `026`: Pandrosion prime atlases, completed zeta, and phase
  destruction programs;
- `027`: Pandrosion relativity as finite-slope causal probes and dynamic
  atlases;
- `028`: numerical Pandrosion relativity benchmarks: Minkowski, Newtonian
  shadow, and horizon-safe atlas switching.

Rebuild a paper with the local Tectonic binary:

```bash
./tectonic latex/tex/028_pandrosion_relativity_benchmarks.tex
```

Copy a rebuilt PDF into the PDF folder when needed:

```bash
cp latex/tex/028_pandrosion_relativity_benchmarks.pdf \
   latex/pdf/028_pandrosion_relativity_benchmarks.pdf
```

### Python engines

The active standalone engines are in [`flow/`](flow/):

- `118` to `120`: probe-aware Thales/Pandrosion engines and equal-value/tensor
  aware NumPy variants;
- `200`, `201`, `215`: backbone, log-stable, and shared core engine work;
- `300` to `304`: atlas, geodesic, basin, universal-atlas, and full-cubic
  Halley NumPy engines;
- `306`: universal atlas plus hypercube inverse-jet correction;
- `307`: full Schwarzschild Pandrosion geodesic solver.

Example: run the 304 autonomous universal-atlas engine:

```bash
.venv/bin/python flow/304_pandrosion_universal_atlas_full_cubic_halley_numpy_engine.py \
  --self-test
```

Example: run the 307 Schwarzschild geodesic solver:

```bash
.venv/bin/python flow/307_full_schwarzschild_pandrosion_geodesic_solver.py \
  --scenario all \
  --out /private/tmp/307_all.json
```

The 307 solver is standalone NumPy code.  It uses Painleve-Gullstrand
coordinates for horizon-safe evolution, dynamic Schwarzschild/PG chart
diagnostics, Hamiltonian-shell projection, finite-slope geodesic defects, and
JSON reporting.

### Figure and benchmark scripts

Recent figure scripts are under [`scripts/`](scripts/):

```bash
.venv/bin/python scripts/025_generate_prime_atlas_figures.py
.venv/bin/python scripts/026_generate_completed_prime_atlas_figures.py
.venv/bin/python scripts/027_generate_pandrosion_relativity_figures.py
.venv/bin/python scripts/028_generate_pandrosion_relativity_benchmarks.py
```

The generated files are written under [`latex/figures/`](latex/figures/).

### Lean formalization

The Lean project lives in [`lean/`](lean/) and uses the pinned toolchain in
[`lean/lean-toolchain`](lean/lean-toolchain).

Run the strict Docker check:

```bash
docker compose run --rm lean-check
```

Run the faster incremental Lean build:

```bash
docker compose run --rm lean-incremental
```

The strict check script audits for `sorry`, `admit`, and raw `axiom`, and
checks theorem dependencies against the local whitelist.

## Repository Layout

```text
.
|-- articles/       # distributable historical PDFs
|-- benchmarks/     # benchmark outputs and local experiment data
|-- flow/           # standalone Python/NumPy Pandrosion engines
|-- latex/
|   |-- figures/    # generated PDF/PNG figures
|   |-- pdf/        # compiled research papers
|   `-- tex/        # LaTeX sources
|-- lean/           # Lean 4 formalization
|-- scripts/        # figure and utility scripts
|-- Dockerfile
|-- docker-compose.yml
|-- LICENSE
|-- README.md
`-- tectonic        # local Tectonic binary
```

## Scope And Non-Claims

This repository deliberately mixes proof work, research exposition, and
experimental code.  They do not have the same evidential status.

- Lean files are compiler-checked formal artifacts.
- LaTeX papers may contain proved statements, conjectures, and research
  programs; check each paper's status paragraph.
- Python scripts are reproducible numerical experiments and solver prototypes.
- Zeta, Smale 17, and physics-facing papers are research programs unless a
  theorem is explicitly proved.

The project does not claim to solve classical open problems or replace standard
physics.  The intended claim is narrower: Pandrosion supplies a finite-slope
geometric calculus with chart scaling, reciprocal/inversion logic, dynamic
atlases, and reproducible numerical experiments.

## Citation

If you use this repository, cite:

```text
Ivan Besevic, Universitas Pandrosion.
DOI: 10.5281/zenodo.19757311
```

DOI landing page: <https://doi.org/10.5281/zenodo.19757311>

## License

This repository is released under the MIT License.  See [`LICENSE`](LICENSE).
