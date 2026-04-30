# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `3,4`
- families: `ks, ks_sparse, dense_iid, sparse_iid, fewnomial, degree_shell_ks, mixed_degree, real_ks, phase_ks, diagonal, chain, cyclic, ill_scaled`
- count/pool/epochs: `4/512/24`
- same generated systems: `True`

## Summary

- 116(115) roots: `48/52`; complete families: `12/13`
- 119(118) roots: `52/52`; complete families: `13/13`
- wall seconds: 116(115)=`0.6317`, 119(118)=`0.4559`, total=`1.0881`
- root delta 119-116: `4`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `3,4` | `ks` | 4 | 4 | 0 | 0.0256 | 0.0353 | 0.727 | 649 | 1895 | `116(115)` |
| `3,4` | `ks_sparse` | 4 | 4 | 0 | 0.0165 | 0.0080 | 2.053 | 692 | 463 | `119(118)` |
| `3,4` | `dense_iid` | 4 | 4 | 0 | 0.0045 | 0.0337 | 0.133 | 198 | 1867 | `116(115)` |
| `3,4` | `sparse_iid` | 4 | 4 | 0 | 0.0143 | 0.0353 | 0.405 | 650 | 1969 | `116(115)` |
| `3,4` | `fewnomial` | 4 | 4 | 0 | 0.0152 | 0.0348 | 0.438 | 704 | 1975 | `116(115)` |
| `3,4` | `degree_shell_ks` | 4 | 4 | 0 | 0.0183 | 0.0595 | 0.307 | 757 | 2201 | `116(115)` |
| `3,4` | `mixed_degree` | 4 | 4 | 0 | 0.0200 | 0.0384 | 0.521 | 890 | 2140 | `116(115)` |
| `3,4` | `real_ks` | 4 | 4 | 0 | 0.0167 | 0.0385 | 0.434 | 750 | 2139 | `116(115)` |
| `3,4` | `phase_ks` | 4 | 4 | 0 | 0.0185 | 0.0349 | 0.531 | 721 | 1919 | `116(115)` |
| `3,4` | `diagonal` | 4 | 4 | 0 | 0.0128 | 0.0184 | 0.693 | 627 | 1076 | `116(115)` |
| `3,4` | `chain` | 4 | 4 | 0 | 0.0157 | 0.0353 | 0.446 | 720 | 2005 | `116(115)` |
| `3,4` | `cyclic` | 4 | 4 | 0 | 0.0146 | 0.0334 | 0.437 | 631 | 1892 | `116(115)` |
| `3,4` | `ill_scaled` | 0 | 4 | 4 | 0.4390 | 0.0503 | 8.722 | 23309 | 2648 | `119(118)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `3,4` | `ks` | 2.135e-09 | 3.665e-13 |
| `3,4` | `ks_sparse` | 3.014e-10 | 9.885e-14 |
| `3,4` | `dense_iid` | 1.635e-09 | 5.636e-12 |
| `3,4` | `sparse_iid` | 5.811e-10 | 2.108e-15 |
| `3,4` | `fewnomial` | 2.990e-09 | 1.421e-15 |
| `3,4` | `degree_shell_ks` | 1.395e-09 | 6.113e-16 |
| `3,4` | `mixed_degree` | 4.420e-10 | 1.143e-11 |
| `3,4` | `real_ks` | 1.880e-09 | 1.049e-14 |
| `3,4` | `phase_ks` | 1.004e-09 | 8.686e-14 |
| `3,4` | `diagonal` | 1.828e-09 | 2.940e-13 |
| `3,4` | `chain` | 1.582e-09 | 1.294e-16 |
| `3,4` | `cyclic` | 3.569e-09 | 4.923e-14 |
| `3,4` | `ill_scaled` |  | 3.840e-15 |
