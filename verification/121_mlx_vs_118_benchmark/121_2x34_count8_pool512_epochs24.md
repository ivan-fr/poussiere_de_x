# 121 MLX vs 118 Benchmark

- cases: `2,34`
- count/pool/epochs: `8/512/24`
- MLX device: `gpu`
- MLX complex dtype: `complex64`

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 121 roots: `0/8`; complete cases: `0/1`
- wall seconds: 118=`0.1816`, 121=`19.9865`, total=`20.1682`
- root delta 121-118: `-8`

## Pair Results

| case | 118 roots | 121 roots | delta | 118 sec | 121 sec | 121 speedup | 118 evals | 121 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 0 | -8 | 0.1816 | 19.9865 | 0.009 | 6535 | 43934 | `118` |

## Residual Minima

| case | 118 min residual | 121 min residual |
|---|---:|---:|
| `2,34` | 2.411e-14 |  |
