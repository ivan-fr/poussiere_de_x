# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `8192/8192/512`
- count/pool/epochs: `32/2048/24`
- scope: Metal batch start selection, Metal n=2 polish pre-pass, 118 complex128 refinement

## Summary

- 118 roots: `32/32`; complete cases: `1/1`
- 123 roots: `32/32`; complete cases: `1/1`
- wall seconds: 118=`2.3056`, 123=`1.3693`, total=`3.6788`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 32 | 32 | 0 | 2.3056 | 1.3693 | 1.684 | 57638 | 17884 | `123(hybrid)` |

## Timing Breakdown

| case | engine | generation | metal select process | metal polish2 process | metal probe process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0137 | 0.0000 | 0.0000 | 0.0000 | 2.2917 | 1.0403 | 0.2415 |
| `2,34` | `123(hybrid)` | 0.0057 | 0.0289 | 0.0718 | 0.0000 | 0.6509 | 0.3057 | 0.0638 |
