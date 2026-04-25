# Numerical Scripts

Run these from the repository root or from `latex/` with the local virtualenv.

Examples:

```bash
.venv/bin/python latex/scripts/verify_monotonicity.py
(cd latex && ../.venv/bin/python scripts/benchmark_general_poly.py)
```

Naming convention:

- `benchmark_*.py`: timing and convergence benchmarks.
- `explore_*.py`: parameter sweeps and exploratory searches.
- `prove_*.py`: symbolic or empirical proof-certificate support.
- `figure*.py` / `figures*.py`: visual asset generators.
- `test_*.py`: small local stress tests.
