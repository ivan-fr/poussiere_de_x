# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `8192/8192/512`
- count/pool/epochs: `32/2048/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `32/32`; complete cases: `1/1`
- 123 roots: `32/32`; complete cases: `1/1`
- wall seconds: 118=`2.5545`, 123=`1.6993`, total=`4.2577`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 32 | 32 | 0 | 2.5545 | 1.6993 | 1.503 | 57638 | 27242 | `123(hybrid)` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0140 | 0.0000 | 2.5403 | 1.1702 | 0.2647 |
| `2,34` | `123(hybrid)` | 0.0025 | 0.0316 | 1.0439 | 0.4944 | 0.1020 |
