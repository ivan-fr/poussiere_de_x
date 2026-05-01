# 128 Swift Metal Port vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `8192/8192/512`
- count/pool/epochs: `32/2048/24`
- scope: same 118 Kostlan system; Metal selector; Metal n=2 polish; Swift Double final Pandrosion refine

## Summary

- 118 roots: `32/32`; complete cases: `1/1`
- 128 roots: `32/32`; complete cases: `1/1`
- wall seconds: 118=`2.3111`, 128=`0.2028`, total=`2.5563`
- root delta 128-118: `0`

## Pair Results

| case | 118 roots | 128 roots | delta | 118 sec | 128 sec | 128 speedup | 118 evals | 128 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 32 | 32 | 0 | 2.3111 | 0.2028 | 11.395 | 57638 | 17097 | `128(swift)` |

## Timing Breakdown

| case | engine | metal select | metal polish2 | refine/extract | eval sec | slope sec | max residual |
|---|---|---:|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0000 | 0.0000 | 2.3107 | 1.0585 | 0.2435 | 6.342e-09 |
| `2,34` | `128(swift)` | 0.0074 | 0.1132 | 0.0489 | 0.0374 | 0.0051 | 7.920e-09 |
