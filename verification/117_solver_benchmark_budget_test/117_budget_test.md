# 117 solver benchmark

- cases: `3,4`
- families: `ks`
- solvers: `scipy_multistart`
- run_external: `True`

## Availability

| solver | available | detail |
|---|---:|---|
| `scipy_multistart` | `True` | SciPy local baseline |

## Results

| case | family | solver | status | roots | evals | seconds | residual min | note |
|---|---|---|---|---:|---:|---:|---:|---|
| `3,4` | `ks` | `scipy_multistart` | `ok` | 0/4 | 10/10 | 0.2492 |  | eval_budget_exhausted |
