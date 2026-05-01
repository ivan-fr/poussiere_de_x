# 123 Swift Metal Hybrid vs 118 Benchmark

- cases: `2,34`
- metal candidates/selector top/refine top: `4096/4096/128`
- count/pool/epochs: `8/512/24`
- scope: Metal batch start selection, 118 complex128 refinement

## Summary

- 118 roots: `8/8`; complete cases: `1/1`
- 123 roots: `8/8`; complete cases: `1/1`
- wall seconds: 118=`0.1992`, 123=`0.3932`, total=`0.5925`
- root delta 123-118: `0`

## Pair Results

| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `2,34` | 8 | 8 | 0 | 0.1992 | 0.3932 | 0.507 | 6535 | 5618 | `118` |

## Timing Breakdown

| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |
|---|---|---:|---:|---:|---:|---:|
| `2,34` | `118` | 0.0095 | 0.0000 | 0.1897 | 0.0856 | 0.0208 |
| `2,34` | `123(hybrid)` | 0.1289 | 0.0550 | 0.1598 | 0.0754 | 0.0158 |
