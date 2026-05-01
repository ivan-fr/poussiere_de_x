# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/refine top: `8192/192`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `3/8`; complete cases: `0/1`
- wall seconds: 118=`0.1842`, 123=`0.6430`, total=`0.8273`
- root delta 123-118: `-5`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 3 | -5 | 0.1842 | 0.6430 | 0.287 | 6535 | 12694 | `118` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0095 | 0.0000 | 0.1746 | 0.0792 | 0.0188 |
| `2,34` | `123(hybrid)` | 0.2207 | 0.0907 | 0.3204 | 0.1534 | 0.0303 |
