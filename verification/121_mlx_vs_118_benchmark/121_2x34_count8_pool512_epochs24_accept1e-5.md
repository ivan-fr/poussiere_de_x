# 121 MLX vs 118 Benchmark

- cases: `2,34`
- count/pool/epochs: `8/512/24`
- MLX device: `gpu`
- MLX complex dtype: `complex64`

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 121 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.0298`, 121=`0.3801`, total=`0.4099`
- root delta 121-118: `0`

## Pair Results

| case | 118 roots | 121 roots | delta | 118 sec | 121 sec | 121 speedup | 118 evals | 121 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.0298 | 0.3801 | 0.078 | 764 | 764 | `118` |

## Residual Minima

| case | 118 min residual | 121 min residual |
|---|---:|---:|
| `2,34` | 3.486e-11 | 8.086e-08 |
