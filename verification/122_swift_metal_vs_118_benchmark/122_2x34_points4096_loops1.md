# 122 Swift Metal vs 118 Benchmark

- cases: `2,34`
- metal points/loops: `4096/1`
- scope: Swift + Metal batch `F(z)` evaluator for the same 118 dense Kostlan system

## Results

| case | terms | 118 eval sec | Metal kernel sec | Metal process sec | kernel speedup | process speedup | full 118 roots | full 118 sec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `2,34` | 630 | 0.0570 | 0.0040 | 0.0518 | 14.124 | 1.101 |  |  |

## Notes

- 122 uses `float2` complex arithmetic in Metal, so it is a batch evaluator prototype, not the final `complex128` corrector.
- The useful number is the process speedup when input/output and command dispatch are included.
- A complete Swift/Metal solver would need to keep probe scoring and candidate selection on the GPU, then refine accepted candidates in CPU `complex128`.
