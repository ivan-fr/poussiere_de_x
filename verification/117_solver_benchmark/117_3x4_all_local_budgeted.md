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

| case | family | solver | status | roots | evals | seconds | residual min | note |
|---|---|---|---|---:|---:|---:|---:|---|
| `3,4` | `ks` | `pandrosion116` | `ok` | 4/4 | 649 | 0.0153 | 2.135e-09 |  |
| `3,4` | `ks` | `scipy_multistart` | `ok` | 3/4 | 512/512 | 0.1422 | 8.066e-11 | eval_budget_exhausted |
| `3,4` | `ks_sparse` | `pandrosion116` | `ok` | 4/4 | 692 | 0.0156 | 3.014e-10 |  |
| `3,4` | `ks_sparse` | `scipy_multistart` | `ok` | 4/4 | 474/512 | 0.0036 | 1.627e-11 | target_count_reached |
| `3,4` | `dense_iid` | `pandrosion116` | `ok` | 4/4 | 198 | 0.0045 | 1.635e-09 |  |
| `3,4` | `dense_iid` | `scipy_multistart` | `ok` | 3/4 | 512/512 | 0.0038 | 9.163e-13 | eval_budget_exhausted |
| `3,4` | `sparse_iid` | `pandrosion116` | `ok` | 4/4 | 650 | 0.0140 | 5.811e-10 |  |
| `3,4` | `sparse_iid` | `scipy_multistart` | `ok` | 2/4 | 512/512 | 0.0040 | 9.515e-11 | eval_budget_exhausted |
| `3,4` | `fewnomial` | `pandrosion116` | `ok` | 4/4 | 704 | 0.0157 | 2.990e-09 |  |
| `3,4` | `fewnomial` | `scipy_multistart` | `ok` | 3/4 | 512/512 | 0.0036 | 2.756e-10 | eval_budget_exhausted |
| `3,4` | `degree_shell_ks` | `pandrosion116` | `ok` | 4/4 | 757 | 0.0166 | 1.395e-09 |  |
| `3,4` | `degree_shell_ks` | `scipy_multistart` | `ok` | 2/4 | 512/512 | 0.0038 | 1.944e-10 | eval_budget_exhausted |
| `3,4` | `mixed_degree` | `pandrosion116` | `ok` | 4/4 | 890 | 0.0206 | 4.420e-10 |  |
| `3,4` | `mixed_degree` | `scipy_multistart` | `ok` | 3/4 | 512/512 | 0.0041 | 3.548e-11 | eval_budget_exhausted |
| `3,4` | `real_ks` | `pandrosion116` | `ok` | 4/4 | 750 | 0.0170 | 1.880e-09 |  |
| `3,4` | `real_ks` | `scipy_multistart` | `ok` | 2/4 | 512/512 | 0.0038 | 1.079e-09 | eval_budget_exhausted |
| `3,4` | `phase_ks` | `pandrosion116` | `ok` | 4/4 | 721 | 0.0166 | 1.004e-09 |  |
| `3,4` | `phase_ks` | `scipy_multistart` | `ok` | 3/4 | 512/512 | 0.0040 | 9.305e-11 | eval_budget_exhausted |
| `3,4` | `diagonal` | `pandrosion116` | `ok` | 4/4 | 627 | 0.0128 | 1.828e-09 |  |
| `3,4` | `diagonal` | `scipy_multistart` | `ok` | 2/4 | 512/512 | 0.0038 | 8.876e-11 | eval_budget_exhausted |
| `3,4` | `chain` | `pandrosion116` | `ok` | 4/4 | 720 | 0.0164 | 1.582e-09 |  |
| `3,4` | `chain` | `scipy_multistart` | `ok` | 4/4 | 384/512 | 0.0029 | 3.236e-10 | target_count_reached |
| `3,4` | `cyclic` | `pandrosion116` | `ok` | 4/4 | 631 | 0.0150 | 3.569e-09 |  |
| `3,4` | `cyclic` | `scipy_multistart` | `ok` | 4/4 | 476/512 | 0.0036 | 3.965e-11 | target_count_reached |
| `3,4` | `ill_scaled` | `pandrosion116` | `ok` | 0/4 | 23309 | 0.4470 |  |  |
| `3,4` | `ill_scaled` | `scipy_multistart` | `ok` | 3/4 | 512/512 | 0.0039 | 3.390e-13 | eval_budget_exhausted |
