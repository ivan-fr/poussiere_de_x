# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `8192/1536/192`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `3/8`; complete cases: `0/1`
- wall seconds: 118=`0.2253`, 123=`0.7152`, total=`0.9405`
- root delta 123-118: `-5`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 3 | -5 | 0.2253 | 0.7152 | 0.315 | 6535 | 12783 | `118` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0318 | 0.0000 | 0.1932 | 0.0862 | 0.0210 |
| `2,34` | `123(hybrid)` | 0.2519 | 0.0685 | 0.3731 | 0.1765 | 0.0354 |
