# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `2048/2048/96`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, Metal n=2 polish pre-pass, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.2449`, 123=`0.1998`, total=`0.4487`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.2449 | 0.1998 | 1.226 | 6535 | 2361 | `123(hybrid)` |

## Timing Breakdown

| case | engine | generation | metal select process | metal polish2 process | metal probe process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0138 | 0.0000 | 0.0000 | 0.0000 | 0.2309 | 0.1034 | 0.0252 |
| `2,34` | `123(hybrid)` | 0.0015 | 0.0078 | 0.0690 | 0.0000 | 0.0890 | 0.0415 | 0.0086 |
