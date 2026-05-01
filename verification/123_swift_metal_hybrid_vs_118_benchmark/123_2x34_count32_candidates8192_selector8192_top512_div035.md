# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `8192/8192/512`
- count/pool/epochs: `32/2048/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `32/32`; complete cases: `1/1`
- 123 roots: `32/32`; complete cases: `1/1`
- wall seconds: 118=`1.7945`, 123=`1.9095`, total=`3.7041`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 32 | 32 | 0 | 1.7945 | 1.9095 | 0.940 | 57638 | 33772 | `118` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0111 | 0.0000 | 1.7833 | 0.8013 | 0.1881 |
| `2,34` | `123(hybrid)` | 0.2633 | 0.0759 | 1.0690 | 0.4958 | 0.1067 |
