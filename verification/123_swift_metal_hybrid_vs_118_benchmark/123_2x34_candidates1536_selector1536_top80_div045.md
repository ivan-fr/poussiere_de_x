# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `1536/1536/80`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.2006`, 123=`0.1981`, total=`0.3988`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.2006 | 0.1981 | 1.013 | 6535 | 3112 | `123(hybrid)` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0109 | 0.0000 | 0.1896 | 0.0858 | 0.0208 |
| `2,34` | `123(hybrid)` | 0.0479 | 0.0420 | 0.0927 | 0.0438 | 0.0093 |
