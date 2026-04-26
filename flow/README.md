# Pandrosion Flow — Python-First Corpus

Each `NNN_<topic>.py` corresponds to one paper in the Pandrosion corpus.
Running the script:

1. Verifies the paper's main numerical claims.
2. Documents the theory / proofs as Python comments and docstrings.

## Numbering convention

- `000` — `pth.pdf` (Pandrosion p-th roots, founding article in `articles/`)
- `001`–`009` — Pandrosion-Smale univariate series (`articles/Npandrosion_smale.pdf`)
- `010` — `9pandrosion_smale_mv.pdf` (multivariate)
- `011`–`100` — Pandrosion legacy series (`articles_pandrosion_legacy/NNpandrosion_*.pdf`)
- `101`–`115+` — Algorithmic + open-conjecture series (`articles_pandrosion_legacy/NNN_*.pdf`)

## Running everything

```bash
cd /Users/ivanbesevic/Documents/poussiere/flow
python3 -m pytest *.py -v   # if pytest-style tests are present
# or
for f in *.py; do echo "=== $f ==="; python3 $f; done
```

## Conventions

Each script header lists:
- `# PAPER: <number>`
- `# TITLE: <title>`
- `# STATUS: <proved | conjectural | empirical | reformulation>`
- `# DEPENDS: <comma-separated list of other paper numbers>`

Each script has at minimum:
- `def verify():` — runs numerical certificate(s).
- `def main():` — orchestrates verification and prints a report.
- `if __name__ == "__main__": main()` at the bottom.

Theory and proofs go in the docstring at the top of the file and inline
near the relevant code.
