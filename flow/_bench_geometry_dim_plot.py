"""Plot geometry-dim sweep results from /tmp/bench_geom_dim_results.json."""
from __future__ import annotations
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

with open("/tmp/bench_geom_dim_results.json") as f:
    results = json.load(f)

# Group by (geom, dim, Tk)
by_key = {}
for r in results:
    key = (r["geom"], r["dim"], r["Tk"])
    by_key.setdefault(key, []).append(r)

# Build per-geom dim curves
GEOM_COLORS = {"path": "tab:blue", "cube": "tab:orange",
               "simplex": "tab:green", "sparse": "tab:red"}
GEOM_LABEL = {"path": "path (T_k order)",
              "cube": "cube (max_orders)",
              "simplex": "simplex (k pivots)",
              "sparse": "sparse (activity power)"}


def avg(rs, k):
    return sum(r[k] for r in rs) / len(rs)


fig, axes = plt.subplots(2, 2, figsize=(13, 9))

for geom in ["path", "cube", "simplex", "sparse"]:
    if geom == "path":
        # path uses T_k as the dimension; key (path, 1, T_k) -> dim is T_k
        keys = sorted([k for k in by_key if k[0] == "path"], key=lambda k: k[2])
        dims = [k[2] for k in keys]
    else:
        keys = sorted([k for k in by_key if k[0] == geom], key=lambda k: k[1])
        dims = [k[1] for k in keys]
    if not keys:
        continue
    cov = [100 * avg(by_key[k], "coverage") for k in keys]
    time_avg = [avg(by_key[k], "time") for k in keys]
    F = [avg(by_key[k], "F_evals") for k in keys]
    coll = [avg(by_key[k], "collisions") for k in keys]
    eff = [(cov[i] / 100.0) / max(F[i], 1) * 1e6 for i in range(len(keys))]

    color = GEOM_COLORS[geom]
    axes[0, 0].plot(dims, cov, "-o", label=GEOM_LABEL[geom], color=color, lw=2)
    axes[0, 1].plot(dims, time_avg, "-o", label=GEOM_LABEL[geom], color=color, lw=2)
    axes[1, 0].plot(dims, F, "-o", label=GEOM_LABEL[geom], color=color, lw=2)
    axes[1, 1].plot(dims, eff, "-o", label=GEOM_LABEL[geom], color=color, lw=2)

axes[0, 0].set_title("Coverage vs dimension (n=3, d=3, Bezout=27)")
axes[0, 0].set_xlabel("dimension"); axes[0, 0].set_ylabel("coverage (%)")
axes[0, 0].grid(alpha=0.3); axes[0, 0].legend()

axes[0, 1].set_title("Wall-clock vs dimension")
axes[0, 1].set_xlabel("dimension"); axes[0, 1].set_ylabel("time (s)")
axes[0, 1].set_yscale("log"); axes[0, 1].grid(alpha=0.3); axes[0, 1].legend()

axes[1, 0].set_title("F-evaluations vs dimension")
axes[1, 0].set_xlabel("dimension"); axes[1, 0].set_ylabel("F-evals")
axes[1, 0].set_yscale("log"); axes[1, 0].grid(alpha=0.3); axes[1, 0].legend()

axes[1, 1].set_title("Efficiency = coverage / F-evals (×10⁶)")
axes[1, 1].set_xlabel("dimension"); axes[1, 1].set_ylabel("efficiency")
axes[1, 1].grid(alpha=0.3); axes[1, 1].legend()

plt.suptitle("Universal-geometry dimension sweep on (n=3,d=3) Kostlan-Smale, "
             "T_k corrector, single-geometry orbit (no portfolio)", fontsize=11)
plt.tight_layout()
out = "/sessions/festive-great-meitner/mnt/poussiere/flow/bench_geometry_dim.png"
plt.savefig(out, dpi=140, bbox_inches="tight")
print(f"saved {out}")
