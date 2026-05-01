# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `2,34`
- families: `ks`
- count/pool/epochs: `3/96/12`
- same generated systems: `True`

## Summary

- 116(115) roots: `1/3`; complete families: `0/1`
- 119(118) roots: `3/3`; complete families: `1/1`
- wall seconds: 116(115)=`0.3612`, 119(118)=`0.0635`, total=`0.4248`
- root delta 119-116: `2`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | `ks` | 1 | 3 | 2 | 0.3612 | 0.0635 | 5.686 | 3710 | 1386 | `119(118)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `2,34` | `ks` | 3.834e-09 | 1.324e-13 |
