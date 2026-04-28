# 076 Pandrosion homothetic geometry package

Run example:

```bash
python -S 076_pandrosion_homothetic_geometry.py \
  --family ks \
  --case 2,8 \
  --retries 2 \
  --base-chunk-size 32 \
  --parallel-base 2 \
  --micro-batch 2 \
  --micro-limit 16 \
  --stop-at-bezout \
  --homothety system \
  --equation-normalize \
  --csv 076_ks28.csv \
  --batch-csv 076_ks28_batches.csv \
  --md 076_ks28.md \
  --roots-json 076_ks28_roots.json \
  --include-lairez-reference
```

The 076 file imports 070, which imports 069, which imports 068. Those dependency files are included here.

## Flow 077: broad benchmark vs Lairez-style

```bash
python -S flow/077_pandrosion_vs_lairez_benchmark.py \
  --suite high \
  --config-set production \
  --seeds 1 \
  --max-bezout 256 \
  --csv 077_benchmark.csv \
  --json 077_benchmark.json \
  --md 077_benchmark.md
```

Useful presets:

- `--suite smoke|quick|high|full` selects KS/dense/sparse systems by Bezout scale.
- `--config-set production|ablation|full` selects the 076 homothety/equation-normalization variants.
- `--family ks --cases 2,8 3,4` overrides the suite with explicit systems.
