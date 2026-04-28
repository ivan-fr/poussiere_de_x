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
