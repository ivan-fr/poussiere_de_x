# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `8192/8192/384`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.1996`, 123=`0.8172`, total=`1.0170`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.1996 | 0.8172 | 0.244 | 6535 | 6598 | `118` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0112 | 0.0000 | 0.1884 | 0.0846 | 0.0204 |
| `2,34` | `123(hybrid)` | 0.2555 | 0.0682 | 0.1950 | 0.0914 | 0.0189 |
