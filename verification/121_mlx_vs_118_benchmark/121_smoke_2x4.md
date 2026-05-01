# 121 MLX vs 118 Benchmark

- cases: `2,4`
- count/pool/epochs: `1/8/2`
- MLX device: `gpu`
- MLX complex dtype: `complex64`

## Summary

- 118 roots: `0/1`; complete cases: `0/1`
- 121 roots: `0/1`; complete cases: `0/1`
- wall seconds: 118=`0.0265`, 121=`0.3058`, total=`0.3324`
- root delta 121-118: `0`

## Pair Results

| case | 118 roots | 121 roots | delta | 118 sec | 121 sec | 121 speedup | 118 evals | 121 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,4` | 0 | 0 | 0 | 0.0265 | 0.3058 | 0.087 | 292 | 292 | `118` |

## Residual Minima

| case | 118 min residual | 121 min residual |
|---|---:|---:|
| `2,4` |  |  |
