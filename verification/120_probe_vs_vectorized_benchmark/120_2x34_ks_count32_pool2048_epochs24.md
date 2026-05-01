# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `2,34`
- families: `ks`
- count/pool/epochs: `32/2048/24`
- same generated systems: `True`

## Summary

- 116(115) roots: `26/32`; complete families: `0/1`
- 119(118) roots: `32/32`; complete families: `1/1`
- wall seconds: 116(115)=`5.7158`, 119(118)=`2.0362`, total=`7.7524`
- root delta 119-116: `6`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | `ks` | 26 | 32 | 6 | 5.7158 | 2.0362 | 2.807 | 102939 | 49184 | `119(118)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `2,34` | `ks` | 1.431e-09 | 1.507e-15 |
