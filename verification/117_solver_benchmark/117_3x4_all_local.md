# 117 solver benchmark

- cases: `3,4`
- families: `ks, ks_sparse, dense_iid, sparse_iid, fewnomial, degree_shell_ks, mixed_degree, real_ks, phase_ks, diagonal, chain, cyclic, ill_scaled`
- solvers: `pandrosion116, scipy_multistart`
- run_external: `True`

## Availability

| solver | available | detail |
|---|---:|---|
| `pandrosion116` | `True` | local Python engine |
| `scipy_multistart` | `True` | SciPy local baseline |

## Results

| case | family | solver | status | roots | seconds | residual min | note |
|---|---|---|---|---:|---:|---:|---|
| `3,4` | `ks` | `pandrosion116` | `ok` | 4/4 | 0.0154 | 2.135e-09 |  |
| `3,4` | `ks` | `scipy_multistart` | `ok` | 4/4 | 0.1343 | 8.066e-11 |  |
| `3,4` | `ks_sparse` | `pandrosion116` | `ok` | 4/4 | 0.0154 | 3.014e-10 |  |
| `3,4` | `ks_sparse` | `scipy_multistart` | `ok` | 4/4 | 0.0035 | 1.627e-11 |  |
| `3,4` | `dense_iid` | `pandrosion116` | `ok` | 4/4 | 0.0045 | 1.635e-09 |  |
| `3,4` | `dense_iid` | `scipy_multistart` | `ok` | 4/4 | 0.0043 | 9.163e-13 |  |
| `3,4` | `sparse_iid` | `pandrosion116` | `ok` | 4/4 | 0.0144 | 5.811e-10 |  |
| `3,4` | `sparse_iid` | `scipy_multistart` | `ok` | 4/4 | 0.0063 | 9.515e-11 |  |
| `3,4` | `fewnomial` | `pandrosion116` | `ok` | 4/4 | 0.0155 | 2.990e-09 |  |
| `3,4` | `fewnomial` | `scipy_multistart` | `ok` | 4/4 | 0.0042 | 2.756e-10 |  |
| `3,4` | `degree_shell_ks` | `pandrosion116` | `ok` | 4/4 | 0.0166 | 1.395e-09 |  |
| `3,4` | `degree_shell_ks` | `scipy_multistart` | `ok` | 4/4 | 0.0055 | 1.880e-10 |  |
| `3,4` | `mixed_degree` | `pandrosion116` | `ok` | 4/4 | 0.0200 | 4.420e-10 |  |
| `3,4` | `mixed_degree` | `scipy_multistart` | `ok` | 4/4 | 0.0069 | 3.548e-11 |  |
| `3,4` | `real_ks` | `pandrosion116` | `ok` | 4/4 | 0.0171 | 1.880e-09 |  |
| `3,4` | `real_ks` | `scipy_multistart` | `ok` | 4/4 | 0.0046 | 5.661e-10 |  |
| `3,4` | `phase_ks` | `pandrosion116` | `ok` | 4/4 | 0.0165 | 1.004e-09 |  |
| `3,4` | `phase_ks` | `scipy_multistart` | `ok` | 4/4 | 0.0044 | 9.305e-11 |  |
| `3,4` | `diagonal` | `pandrosion116` | `ok` | 4/4 | 0.0130 | 1.828e-09 |  |
| `3,4` | `diagonal` | `scipy_multistart` | `ok` | 4/4 | 0.0048 | 8.876e-11 |  |
| `3,4` | `chain` | `pandrosion116` | `ok` | 4/4 | 0.0163 | 1.582e-09 |  |
| `3,4` | `chain` | `scipy_multistart` | `ok` | 4/4 | 0.0028 | 3.236e-10 |  |
| `3,4` | `cyclic` | `pandrosion116` | `ok` | 4/4 | 0.0149 | 3.569e-09 |  |
| `3,4` | `cyclic` | `scipy_multistart` | `ok` | 4/4 | 0.0033 | 3.965e-11 |  |
| `3,4` | `ill_scaled` | `pandrosion116` | `ok` | 0/4 | 0.4425 |  |  |
| `3,4` | `ill_scaled` | `scipy_multistart` | `ok` | 4/4 | 0.0048 | 3.390e-13 |  |
