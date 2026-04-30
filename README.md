# Universitas Pandrosion

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19757311.svg)](https://doi.org/10.5281/zenodo.19757311)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)
[![Universitas Pandrosion CI](https://github.com/ivan-fr/poussiere_de_x/actions/workflows/ci.yml/badge.svg)](https://github.com/ivan-fr/poussiere_de_x/actions)

Pandrosion is a research repository around derivative-free root-finding,
Steffensen acceleration, and multivariate polynomial systems.  The repository
has two distinct parts:

- a Lean 4 formalization under [`lean/`](lean/) for the core Pandrosion theory;
- experimental Python/LaTeX work under [`flow/`](flow/), [`latex/`](latex/),
  [`articles/`](articles/), and [`figures/`](figures/).

The Lean code is the high-rigour part of the repository.  The Python and LaTeX
material is reproducible research code and exposition, not formal proof unless
explicitly stated.

## Current Checkout

This checkout currently contains:

- **121 Lean modules** under [`lean/Pandrosion/`](lean/Pandrosion/), using
  `leanprover/lean4:v4.7.0` and `mathlib4 v4.7.0`;
- **compiled and source LaTeX papers** for the Pth/Smale/Pandrosion line,
  including the current papers on engines 113 and 115;
- **7 active `flow/` Python scripts**, focused on the recent autonomous
  multivariate Pandrosion engines:
  - [`112_pandrosion_multifamily_heuristic_thales_engine.py`](flow/112_pandrosion_multifamily_heuristic_thales_engine.py)
  - [`113_pandrosion_heuristic_thales_engine.py`](flow/113_pandrosion_heuristic_thales_engine.py)
  - [`113_pandrosion_autonomous_heuristic_thales_engine.py`](flow/113_pandrosion_autonomous_heuristic_thales_engine.py)
  - [`115_pandrosion_vectorized_pure_pandrosion.py`](flow/115_pandrosion_vectorized_pure_pandrosion.py)
  - [`115_pandrosion_vectorized_pure_thales_engine.py`](flow/115_pandrosion_vectorized_pure_thales_engine.py)
  - [`116_pandrosion_multifamily_vectorized_pure_pandrosion.py`](flow/116_pandrosion_multifamily_vectorized_pure_pandrosion.py)
  - [`117_pandrosion_solver_benchmark.py`](flow/117_pandrosion_solver_benchmark.py)
- generated figures for engines 115 and 116 in [`figures/`](figures/).

Older README versions described a much larger historical `flow/` corpus.  This
README describes the files present in this checkout.

## Research Status

### Formal core

The Lean project is in [`lean/`](lean/).  The strict check script performs a
clean build, scans for `sorry`, `admit`, and raw `axiom`, and audits theorem
dependencies against the classical whitelist
`{propext, Classical.choice, Quot.sound}`.

Run the authoritative check with Docker:

```bash
docker compose run --rm lean-check
```

For faster iteration, use the incremental build:

```bash
docker compose run --rm lean-incremental
```

The Lean package can also be entered interactively:

```bash
docker compose run --rm lean
```

### Experimental engines

The recent Python engines study direct multistart extraction of complex roots
of square polynomial systems.

- **113** uses the Thales/Riemann/Mobius geometric start layer with a local
  heuristic corrector.
- **115** replaces the local Jacobian-style corrector with an exact finite
  Pandrosion slope matrix `Q(a,b)` built by a monomial telescopic identity and
  vectorized as a dense `C @ T` product.
- **116** reuses the 115 flow and corrector, but expands the input generators
  to many multivariate families: dense Kostlan, sparse Kostlan, IID dense and
  sparse, fewnomial, degree-shell, mixed-degree, real, phase-only, structured
  diagonal/chain/cyclic, and ill-scaled systems.
- **117** is a benchmark harness.  It runs `pandrosion116`, runs a budgeted
  local `scipy.optimize.root` multistart baseline, and exports the same systems
  to Bertini, PHCpack, and Julia HomotopyContinuation.jl formats.  External
  solvers are only run when their commands are installed; otherwise they are
  reported as `skipped`.

The current local results show that 116 performs well on the tested
multifamily benchmarks and beats the budgeted SciPy multistart baseline on
coverage for `ks(3,4) --families all`.  This is **not** yet a claim of
superiority over Bertini, PHCpack, HomotopyContinuation.jl, or Lairez-style
homotopy algorithms.  Script 117 exists to make those comparisons explicit.

## Repository Layout

```text
.
├── articles/       # distributable PDFs
├── figures/        # generated figures for papers and engine summaries
├── flow/           # current Python experiment and benchmark scripts
├── latex/          # LaTeX sources and compiled local PDFs
├── lean/           # Lean 4 formalization
├── scripts/        # utility scripts for figures and Lean iteration
├── Dockerfile
├── docker-compose.yml
├── README.md
└── tectonic        # local Tectonic binary for LaTeX builds
```

## Build The Papers

The repository includes a local `tectonic` binary.  To rebuild the current
paper on engine 115:

```bash
./tectonic latex/paper_115_pandrosion_pure_vectorized.tex
```

To copy the rebuilt PDF into `articles/`:

```bash
cp latex/paper_115_pandrosion_pure_vectorized.pdf \
   articles/paper_115_pandrosion_pure_vectorized.pdf
```

The current engine-115 paper is:

- source: [`latex/paper_115_pandrosion_pure_vectorized.tex`](latex/paper_115_pandrosion_pure_vectorized.tex)
- PDF: [`articles/paper_115_pandrosion_pure_vectorized.pdf`](articles/paper_115_pandrosion_pure_vectorized.pdf)

## Run The Recent Experiments

Use the repository virtual environment when available:

```bash
.venv/bin/python --version
```

### Engine 113

```bash
.venv/bin/python flow/113_pandrosion_heuristic_thales_engine.py \
  --cases 2,10 \
  --outdir verification/113_heuristic_out
```

### Engine 115

```bash
.venv/bin/python flow/115_pandrosion_vectorized_pure_pandrosion.py \
  --cases 2,10 \
  --outdir verification/115_vectorized_pure_out
```

### Engine 116: multifamily generator layer

List available families:

```bash
.venv/bin/python flow/116_pandrosion_multifamily_vectorized_pure_pandrosion.py \
  --list-families
```

Run all families on a multivariate case:

```bash
.venv/bin/python flow/116_pandrosion_multifamily_vectorized_pure_pandrosion.py \
  --cases 3,4 \
  --families all \
  --count 4 \
  --pool 512 \
  --outdir verification/116_multifamily_out
```

### Engine 117: benchmark harness

List solver adapters:

```bash
.venv/bin/python flow/117_pandrosion_solver_benchmark.py --list-solvers
```

Run the local comparison between `pandrosion116` and budgeted SciPy:

```bash
.venv/bin/python flow/117_pandrosion_solver_benchmark.py \
  --cases 3,4 \
  --families all \
  --count 4 \
  --pool 512 \
  --solvers local \
  --scipy-starts 128 \
  --outdir verification/117_solver_benchmark
```

The SciPy baseline is budgeted.  By default, `--scipy-eval-budget 0` means
“use `--pool` as the global number of residual evaluations per case/family.”
Use `--scipy-eval-budget -1` only when you intentionally want an unbounded
baseline.

Run all adapters, including external solver exports:

```bash
.venv/bin/python flow/117_pandrosion_solver_benchmark.py \
  --cases 3,4 \
  --families all \
  --count 4 \
  --pool 512 \
  --solvers all \
  --outdir verification/117_solver_benchmark
```

If `bertini`, `phc`, or `julia` are not in `PATH`, script 117 still writes the
input files and marks those solvers as `skipped`.

## Generate Figures

The figure scripts use Matplotlib and the local virtual environment.

```bash
.venv/bin/python scripts/make_115_pandrosion_figure.py
.venv/bin/python scripts/make_116_multifamily_multivariate_figure.py
```

Generated outputs:

- [`figures/115_pandrosion_vectorized_pure_pandrosion.png`](figures/115_pandrosion_vectorized_pure_pandrosion.png)
- [`figures/115_pandrosion_vectorized_pure_pandrosion.pdf`](figures/115_pandrosion_vectorized_pure_pandrosion.pdf)
- [`figures/116_multifamily_multivariate_all.png`](figures/116_multifamily_multivariate_all.png)
- [`figures/116_multifamily_multivariate_all.pdf`](figures/116_multifamily_multivariate_all.pdf)

## Scope And Non-Claims

This repository contains several kinds of evidence:

- Lean theorem proofs, checked by the Lean compiler and strict audit script;
- mathematical exposition in LaTeX;
- numerical experiments in Python;
- generated figures and benchmark reports.

Only the Lean artifacts are formal proofs.  The Python scripts are
reproducible numerical experiments and benchmark harnesses.  The repository
does not claim to resolve classical open problems unless a statement is
explicitly formalized or otherwise proved in the relevant paper.

For the recent multivariate engines specifically:

- `115` establishes and implements the vectorized exact telescopic slope
  corrector;
- `116` shows that the method is robust across many tested polynomial
  families, with `ill_scaled` deliberately acting as a stress test;
- `117` prepares honest comparisons against standard solvers, but external
  homotopy solvers must actually be installed and run before making any
  publication-level performance claim against them.

## Citation

If you use this repository, cite:

```text
Ivan Besevic, Universitas Pandrosion.
DOI: 10.5281/zenodo.19757311
```

DOI landing page: <https://doi.org/10.5281/zenodo.19757311>

## License

See [`LICENSE`](LICENSE).  The repository badge advertises CC BY-SA 4.0; check
the license file for the authoritative terms included with this checkout.
