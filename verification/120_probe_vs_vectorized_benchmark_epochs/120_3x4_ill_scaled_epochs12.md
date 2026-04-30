# 120 Pandrosion Probe-Aware vs Vectorized Benchmark

- cases: `3,4`
- families: `ill_scaled`
- count/pool/epochs: `4/512/12`
- same generated systems: `True`

## Summary

- 116(115) roots: `0/4`; complete families: `0/1`
- 119(118) roots: `4/4`; complete families: `1/1`
- wall seconds: 116(115)=`0.3879`, 119(118)=`0.0828`, total=`0.4708`
- root delta 119-116: `4`

## Pair Results

| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `3,4` | `ill_scaled` | 0 | 4 | 4 | 0.3879 | 0.0828 | 4.687 | 19437 | 4479 | `119(118)` |

## Residual Minima

| case | family | 116 min residual | 119 min residual |
|---|---|---:|---:|
| `3,4` | `ill_scaled` |  | 6.187e-13 |
