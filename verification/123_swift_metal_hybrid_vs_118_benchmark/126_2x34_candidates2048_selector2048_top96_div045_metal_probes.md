# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `2048/2048/96`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.2102`, 123=`0.6103`, total=`0.8244`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.2102 | 0.6103 | 0.344 | 6535 | 1963 | `118` |

## Timing Breakdown

| case | engine | generation | metal select process | metal probe process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0113 | 0.0000 | 0.0000 | 0.1988 | 0.0899 | 0.0216 |
| `2,34` | `123(hybrid)` | 0.0049 | 0.0066 | 0.4580 | 0.5783 | 0.0281 | 0.0177 |
