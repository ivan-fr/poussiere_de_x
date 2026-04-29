# 117 solver benchmark

- cases: `3,4`
- families: `ks`
- solvers: `pandrosion116, scipy_multistart`
- run_external: `True`

## Availability

| solver | available | detail |
|---|---:|---|
| `pandrosion116` | `True` | local Python engine |
| `scipy_multistart` | `True` | SciPy local baseline |

## Results

| case | family | solver | status | roots | evals | seconds | residual min | note |
|---|---|---|---|---:|---:|---:|---:|---|
| `3,4` | `ks` | `pandrosion116` | `ok` | 2/2 | 531 | 0.0128 | 4.078e-10 |  |
| `3,4` | `ks` | `scipy_multistart` | `ok` | 1/2 | 64/64 | 0.2434 | 8.066e-11 | eval_budget_exhausted |
