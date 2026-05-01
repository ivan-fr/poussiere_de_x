# 128 Swift Metal Port vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `2048/2048/96`
- count/pool/epochs: `8/512/24`
- scope: same 118 Kostlan system; Metal selector; Metal n=2 polish; Swift Double final Pandrosion refine

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 128 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.2343`, 128=`0.0913`, total=`0.3672`
- root delta 128-118: `0`

## Pair Results

| case | 118 roots | 128 roots | delta | 118 sec | 128 sec | 128 speedup | 118 evals | 128 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.2343 | 0.0913 | 2.567 | 6535 | 4205 | `128(swift)` |

## Timing Breakdown

| case | engine | metal select | metal polish2 | refine/extract | eval sec | slope sec | max residual |
|---|---|---:|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0000 | 0.0000 | 0.2340 | 0.1067 | 0.0259 | 1.213e-09 |
| `2,34` | `128(swift)` | 0.0039 | 0.0668 | 0.0113 | 0.0086 | 0.0012 | 3.577e-09 |
