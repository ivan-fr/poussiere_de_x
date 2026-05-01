# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `2048/2048/96`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.2179`, 123=`0.2084`, total=`0.4264`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.2179 | 0.2084 | 1.045 | 6535 | 2485 | `123(hybrid)` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0112 | 0.0000 | 0.2066 | 0.0936 | 0.0227 |
| `2,34` | `123(hybrid)` | 0.0686 | 0.0450 | 0.0715 | 0.0335 | 0.0069 |
