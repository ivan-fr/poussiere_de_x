# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `3,4`
- families: `ill_scaled`
- count/pool/epochs: `4/512/8`
- same generated systems: `True`

## Summary

- 116(115) roots: `0/4`; complete families: `0/1`
- 119(118) roots: `4/4`; complete families: `1/1`
- wall seconds: 116(115)=`0.3367`, 119(118)=`0.0329`, total=`0.3697`
- root delta 119-116: `4`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `3,4` | `ill_scaled` | 0 | 4 | 4 | 0.3367 | 0.0329 | 10.225 | 17170 | 1738 | `119(118)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `3,4` | `ill_scaled` |  | 7.819e-14 |
