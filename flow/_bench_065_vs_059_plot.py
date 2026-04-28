"""Generate plot from bench results in /tmp/bench_065_results.json."""
from __future__ import annotations
import json, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

with open("/tmp/bench_065_results.json") as f:
    results = json.load(f)

# Group by (n, d)
by_nd = {}
for r in results:
    by_nd.setdefault((r["n"], r["d"]), []).append(r)

cases = sorted(by_nd.keys())
labels = [f"({n},{d})\nBez={by_nd[(n,d)][0]['bezout']}" for (n, d) in cases]

cov_065 = [100 * sum(r["cov_065"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]
cov_059 = [100 * sum(r["cov_059"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]
t_065 = [sum(r["time_065"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]
t_059 = [sum(r["time_059"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]
F_065 = [sum(r["F_065"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]
F_059 = [sum(r["F_059"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]
coll_065 = [sum(r["coll_065"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]
coll_059 = [sum(r["coll_059"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]
ret_065 = [sum(r["retries_065"] for r in by_nd[k]) / len(by_nd[k]) for k in cases]

fig, axes = plt.subplots(2, 2, figsize=(13, 9))
x = np.arange(len(cases)); w = 0.4

# Coverage
ax = axes[0, 0]
ax.bar(x - w/2, cov_065, w, label="065 γ-retry", color="steelblue")
ax.bar(x + w/2, cov_059, w, label="059 γ Lairez", color="darkorange")
ax.axhline(100, color="gray", ls="--", alpha=0.5)
ax.set_ylabel("Coverage (%)"); ax.set_title("Coverage Bezout (avg over seeds)")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylim(85, 102); ax.legend(); ax.grid(alpha=0.3)

# Wall-clock time (log scale)
ax = axes[0, 1]
ax.bar(x - w/2, t_065, w, label="065", color="steelblue")
ax.bar(x + w/2, t_059, w, label="059", color="darkorange")
ax.set_ylabel("Wall-clock time (s)"); ax.set_title("Wall-clock per system (log)")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.legend(); ax.grid(alpha=0.3)

# F-evaluations (log)
ax = axes[1, 0]
ax.bar(x - w/2, F_065, w, label="065", color="steelblue")
ax.bar(x + w/2, F_059, w, label="059", color="darkorange")
ax.set_ylabel("F evaluations"); ax.set_title("Polynomial evaluations per system (log)")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.legend(); ax.grid(alpha=0.3)

# Collisions and retries
ax = axes[1, 1]
ax.bar(x - w/2, coll_065, w, label="065 collisions", color="steelblue", alpha=0.7)
ax.bar(x - w/2, ret_065, w, bottom=coll_065, label="065 retries used",
       color="navy", alpha=0.9)
ax.bar(x + w/2, coll_059, w, label="059 collisions", color="darkorange", alpha=0.7)
ax.set_ylabel("Avg paths / retries"); ax.set_title("Collisions and γ-retries")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.legend(fontsize=8); ax.grid(alpha=0.3)

plt.suptitle("Benchmark 065 (γ-retry) vs 059 (Lairez γ fixed) on Kostlan-Smale systems",
             fontsize=12)
plt.tight_layout()
out = "/sessions/festive-great-meitner/mnt/poussiere/flow/bench_065_vs_059.png"
plt.savefig(out, dpi=140, bbox_inches="tight")
print(f"saved {out}")
print(f"\n# instances by (n,d):", {k: len(by_nd[k]) for k in cases})
