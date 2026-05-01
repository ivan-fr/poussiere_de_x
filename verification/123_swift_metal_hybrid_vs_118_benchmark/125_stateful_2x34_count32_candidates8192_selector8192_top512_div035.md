# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `8192/8192/512`
- count/pool/epochs: `32/2048/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `32/32`; complete cases: `1/1`
- 123 roots: `32/32`; complete cases: `1/1`
- wall seconds: 118=`2.2667`, 123=`1.6839`, total=`3.9545`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 32 | 32 | 0 | 2.2667 | 1.6839 | 1.346 | 57638 | 27242 | `123(hybrid)` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0147 | 0.0000 | 2.2518 | 1.0246 | 0.2387 |
| `2,34` | `123(hybrid)` | 0.0038 | 0.0289 | 1.0555 | 0.4977 | 0.1049 |
