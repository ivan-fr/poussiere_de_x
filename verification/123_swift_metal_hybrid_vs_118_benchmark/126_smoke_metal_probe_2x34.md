# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `256/256/16`
- count/pool/epochs: `1/32/4`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `0/1`; complete cases: `0/1`
- 123 roots: `1/1`; complete cases: `1/1`
- wall seconds: 118=`0.0775`, 123=`0.0236`, total=`0.1050`
- root delta 123-118: `1`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 0 | 1 | 1 | 0.0775 | 0.0236 | 3.291 | 1686 | 24 | `123(hybrid)` |

## Timing Breakdown

| case | engine | generation | metal select process | metal probe process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0247 | 0.0000 | 0.0000 | 0.0526 | 0.0242 | 0.0046 |
| `2,34` | `123(hybrid)` | 0.0011 | 0.0028 | 0.0143 | 0.0159 | 0.0004 | 0.0002 |
