# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `8192/8192/512`
- count/pool/epochs: `32/2048/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `32/32`; complete cases: `1/1`
- 123 roots: `32/32`; complete cases: `1/1`
- wall seconds: 118=`2.3052`, 123=`3.3516`, total=`5.6604`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 32 | 32 | 0 | 2.3052 | 3.3516 | 0.688 | 57638 | 10410 | `118` |

## Timing Breakdown

| case | engine | generation | metal select process | metal probe process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0129 | 0.0000 | 0.0000 | 2.2911 | 1.0384 | 0.2395 |
| `2,34` | `123(hybrid)` | 0.0036 | 0.0286 | 2.0435 | 2.7948 | 0.1799 | 0.1103 |
