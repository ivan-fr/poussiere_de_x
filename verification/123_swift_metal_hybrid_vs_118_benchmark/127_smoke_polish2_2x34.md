# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `256/256/16`
- count/pool/epochs: `1/32/4`
- scope: Metal batch start selection, Metal n=2 polish pre-pass, 118 complex128 refinement

## Summary

- 118 roots: `0/1`; complete cases: `0/1`
- 123 roots: `1/1`; complete cases: `1/1`
- wall seconds: 118=`0.0849`, 123=`0.0514`, total=`0.1402`
- root delta 123-118: `1`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 0 | 1 | 1 | 0.0849 | 0.0514 | 1.653 | 1686 | 47 | `123(hybrid)` |

## Timing Breakdown

| case | engine | generation | metal select process | metal polish2 process | metal probe process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0236 | 0.0000 | 0.0000 | 0.0000 | 0.0612 | 0.0284 | 0.0051 |
| `2,34` | `123(hybrid)` | 0.0047 | 0.0030 | 0.0359 | 0.0000 | 0.0018 | 0.0008 | 0.0001 |
