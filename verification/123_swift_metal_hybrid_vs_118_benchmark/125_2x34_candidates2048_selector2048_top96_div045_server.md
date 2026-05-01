# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `2048/2048/96`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.3129`, 123=`0.3320`, total=`0.6488`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.3129 | 0.3320 | 0.942 | 6535 | 5123 | `118` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0243 | 0.0000 | 0.2884 | 0.1286 | 0.0309 |
| `2,34` | `123(hybrid)` | 0.0024 | 0.0874 | 0.2121 | 0.0993 | 0.0207 |
