# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `3,4`
- families: `ill_scaled`
- count/pool/epochs: `4/512/4`
- same generated systems: `True`

## Summary

- 116(115) roots: `0/4`; complete families: `0/1`
- 119(118) roots: `0/4`; complete families: `0/1`
- wall seconds: 116(115)=`0.2629`, 119(118)=`0.5363`, total=`0.7993`
- root delta 119-116: `0`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `3,4` | `ill_scaled` | 0 | 0 | 0 | 0.2629 | 0.5363 | 0.490 | 13657 | 29586 | `116(115)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `3,4` | `ill_scaled` |  |  |
