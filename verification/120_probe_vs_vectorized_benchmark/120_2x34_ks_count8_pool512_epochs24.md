# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `2,34`
- families: `ks`
- count/pool/epochs: `8/512/24`
- same generated systems: `True`

## Summary

- 116(115) roots: `8/8`; complete families: `1/1`
- 119(118) roots: `8/8`; complete families: `1/1`
- wall seconds: 116(115)=`0.6025`, 119(118)=`0.6808`, total=`1.2836`
- root delta 119-116: `0`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | `ks` | 8 | 8 | 0 | 0.6025 | 0.6808 | 0.885 | 6633 | 8526 | `116(115)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `2,34` | `ks` | 2.165e-09 | 6.805e-15 |
