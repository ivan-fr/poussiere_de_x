# 117 solver benchmark

- cases: `3,4`
- families: `ks, diagonal`
- solvers: `pandrosion116, scipy_multistart, bertini, phcpack, julia_hc, lairez_custom`
- run_external: `True`

## Availability

| solver | available | detail |
|---|---:|---|
| `pandrosion116` | `True` | local Python engine |
| `scipy_multistart` | `True` | SciPy local baseline |
| `bertini` | `False` | bertini |
| `phcpack` | `False` | phc |
| `julia_hc` | `False` | julia |
| `lairez_custom` | `False` |  |

## Results

| case | family | solver | status | roots | seconds | residual min | note |
|---|---|---|---|---:|---:|---:|---|
| `3,4` | `ks` | `pandrosion116` | `ok` | 2/2 | 0.0189 | 4.078e-10 |  |
| `3,4` | `ks` | `scipy_multistart` | `ok` | 2/2 | 0.4473 | 8.066e-11 |  |
| `3,4` | `ks` | `bertini` | `skipped` |  | 0.0000 |  | command not found: bertini |
| `3,4` | `ks` | `phcpack` | `skipped` |  | 0.0000 |  | command not found: phc |
| `3,4` | `ks` | `julia_hc` | `skipped` |  | 0.0000 |  | command not found: julia |
| `3,4` | `ks` | `lairez_custom` | `skipped` |  | 0.0000 |  | no standard Lairez CLI is configured; pass --lairez-command with placeholders such as {system_json} |
| `3,4` | `diagonal` | `pandrosion116` | `ok` | 2/2 | 0.0120 | 1.828e-09 |  |
| `3,4` | `diagonal` | `scipy_multistart` | `ok` | 2/2 | 0.0021 | 8.876e-11 |  |
| `3,4` | `diagonal` | `bertini` | `skipped` |  | 0.0000 |  | command not found: bertini |
| `3,4` | `diagonal` | `phcpack` | `skipped` |  | 0.0000 |  | command not found: phc |
| `3,4` | `diagonal` | `julia_hc` | `skipped` |  | 0.0000 |  | command not found: julia |
| `3,4` | `diagonal` | `lairez_custom` | `skipped` |  | 0.0000 |  | no standard Lairez CLI is configured; pass --lairez-command with placeholders such as {system_json} |
