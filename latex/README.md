# LaTeX Workspace

This directory is intentionally split by artifact type:

- `*.tex`: source papers and the combined-edition driver.
- `*.pdf`: compiled paper PDFs used by `universitas_pandrosion_combined.tex`.
- `scripts/`: Python experiments, benchmarks, proof-certificate scripts, and figure generators.
- `figures/`: generated visual assets used by the papers.
- `figs/`: Paper-0 figure set used through `\graphicspath`.
- `archive/`: legacy PDFs not used by the current build.

Build from this directory with:

```bash
../tectonic 0pandrosion_pth.tex
../tectonic universitas_pandrosion_combined.tex
```

Run experiments from `scripts/`, for example:

```bash
../.venv/bin/python scripts/verify_monotonicity.py
```
