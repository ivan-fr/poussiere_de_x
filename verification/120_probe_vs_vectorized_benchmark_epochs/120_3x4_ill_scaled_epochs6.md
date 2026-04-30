# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `3,4`
- families: `ill_scaled`
- count/pool/epochs: `4/512/6`
- same generated systems: `True`

## Summary

- 116(115) roots: `0/4`; complete families: `0/1`
- 119(118) roots: `4/4`; complete families: `1/1`
- wall seconds: 116(115)=`0.2823`, 119(118)=`0.1996`, total=`0.4821`
- root delta 119-116: `4`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `3,4` | `ill_scaled` | 0 | 4 | 4 | 0.2823 | 0.1996 | 1.414 | 15975 | 11441 | `119(118)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `3,4` | `ill_scaled` |  | 2.059e-12 |
