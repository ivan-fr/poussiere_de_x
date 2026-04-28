"""
PAPER: 048
TITLE: Dimension scaling for latex/1 universal hypercube and simplex geometries
STATUS: benchmark requested after flows 046/047

QUESTION
========

What happens when the universal geometries from latex/1pandrosion_smale.tex
are pushed to higher dimensions?

Test:
    - hypercube universal geometry from flow/046 at n=5,6
    - simplex universal geometry from flow/047 at n=4,5,6

All systems use the corresponding closed-form Pandrosion geometry as ground
truth:
    hypercube: P_i=x_i(1-s_i)prod_j S_p(s_j)-(x_i-1)
    simplex:   P_i=x_i(1-s_i)S_simplex(s)-(x_i-1)

Compare against the structural quotient prism from flow/025/027.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_flow(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


U_HYPER = load_flow("flow046_universal_hyper", "046_latex1_universal_hypercube.py")
U_SIMPLEX = load_flow("flow047_universal_simplex", "047_latex1_universal_simplex.py")


def primes(n):
    pool = [2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0]
    return pool[:n]


def tag_structural(result):
    return f"q{result['qdim']} {result['iters']}it/{result['evals']} r={result['res']:.0e}"


def tag_universal(result):
    return f"{result['steps']}st/{result['evals']} r={result['res']:.0e} {'OK' if result['ok'] else 'NO'}"


def tag_multi(result):
    return f"{result['success']}/{result['total']} roots={len(result['roots'])} best={result['best']['evals']}"


def run_hypercube(n, p=2, multistart_count=10):
    x_vec = primes(n)
    start = [0.5] * n
    structural = U_HYPER.structural_025(x_vec, p)
    one_path = U_HYPER.universal_reanchor(start, x_vec, p, mode="path")
    avg_cube = U_HYPER.universal_reanchor(start, x_vec, p, mode="avg")
    multi = U_HYPER.multistart(x_vec, p, mode="avg", count=multistart_count)
    return {
        "n": n,
        "p": p,
        "x_vec": x_vec,
        "structural": structural,
        "one_path": one_path,
        "avg_cube": avg_cube,
        "multi": multi,
    }


def run_simplex(n, p=2, multistart_count=10):
    x_vec = primes(n)
    start = [0.5] * n
    structural = U_SIMPLEX.structural_027(x_vec, p)
    universal = U_SIMPLEX.universal_reanchor(start, x_vec, p)
    multi = U_SIMPLEX.multistart(x_vec, p, count=multistart_count)
    return {
        "n": n,
        "p": p,
        "x_vec": x_vec,
        "structural": structural,
        "universal": universal,
        "multi": multi,
    }


def main():
    print("=" * 118)
    print("flow/048 -- dimension scaling: universal hypercube vs universal simplex")
    print("=" * 118)
    print("n = number of dynamic variables; geometric Pandrosion dimension is n+1.")
    print("x_vec = first n primes; p=2 unless stated.")
    print()

    print("HYPERCUBE universal geometry from flow/046")
    print("-" * 118)
    print(f"{'n':>3} {'geom dim':>8} {'p':>2} | {'025 structural':>22} | {'one path':>24} | {'cube avg':>24} | {'multi avg':>18}")
    for n in [5, 6]:
        row = run_hypercube(n, p=2, multistart_count=10)
        print(
            f"{n:>3} {n+1:>8} {row['p']:>2} | "
            f"{tag_structural(row['structural']):>22} | "
            f"{tag_universal(row['one_path']):>24} | "
            f"{tag_universal(row['avg_cube']):>24} | "
            f"{tag_multi(row['multi']):>18}"
        )

    print()
    print("SIMPLEX universal geometry from flow/047")
    print("-" * 118)
    print(f"{'n':>3} {'geom dim':>8} {'p':>2} | {'027 structural':>22} | {'simplex universal':>24} | {'multistart':>18}")
    for n in [4, 5, 6]:
        row = run_simplex(n, p=2, multistart_count=12)
        print(
            f"{n:>3} {n+1:>8} {row['p']:>2} | "
            f"{tag_structural(row['structural']):>22} | "
            f"{tag_universal(row['universal']):>24} | "
            f"{tag_multi(row['multi']):>18}"
        )

    print()
    print("Extra: simplex p=3 at n=4,5")
    print("-" * 118)
    print(f"{'n':>3} {'geom dim':>8} {'p':>2} | {'027 structural':>22} | {'simplex universal':>24} | {'multistart':>18}")
    for n in [4, 5]:
        row = run_simplex(n, p=3, multistart_count=8)
        print(
            f"{n:>3} {n+1:>8} {row['p']:>2} | "
            f"{tag_structural(row['structural']):>22} | "
            f"{tag_universal(row['universal']):>24} | "
            f"{tag_multi(row['multi']):>18}"
        )

    print()
    print("=" * 118)
    print("READING")
    print("=" * 118)
    print("  - Hypercube avg is the most geometrically faithful latex/1 cube, but")
    print("    averaging many Schmidt paths is expensive as n grows.")
    print("  - One-path hypercube is the practical Schmidt version: order-dependent")
    print("    but much cheaper and still exact.")
    print("  - Simplex universal scales better here because Q_simplex uses n simplex")
    print("    edges plus one radial correction, not many cube paths.")
    print("  - Structural prisms stay cheaper whenever the closed-form F_iter exists.")
    print("=" * 118)


if __name__ == "__main__":
    main()
