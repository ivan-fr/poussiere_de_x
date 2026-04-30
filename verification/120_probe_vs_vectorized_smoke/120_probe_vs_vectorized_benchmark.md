# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `2,3`
- families: `ks, diagonal`
- count/pool/epochs: `2/64/12`
- same generated systems: `True`

## Summary

- 116(115) roots: `4/4`; complete families: `2/2`
- 119(118) roots: `4/4`; complete families: `2/2`
- wall seconds: 116(115)=`0.0322`, 119(118)=`0.0377`, total=`0.0700`
- root delta 119-116: `0`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,3` | `ks` | 2 | 2 | 0 | 0.0304 | 0.0188 | 1.621 | 402 | 1203 | `119(118)` |
| `2,3` | `diagonal` | 2 | 2 | 0 | 0.0018 | 0.0189 | 0.097 | 90 | 1235 | `116(115)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `2,3` | `ks` | 1.988e-09 | 2.575e-10 |
| `2,3` | `diagonal` | 8.949e-10 | 8.460e-14 |
