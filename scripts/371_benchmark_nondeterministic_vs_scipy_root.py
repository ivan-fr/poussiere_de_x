"""Blind paired benchmark: Pandrosion 371 vs scipy.optimize.root(hybr).

The manifest is generated from fresh OS entropy unless ``--seed`` is passed,
then sealed and hashed before any solver run.  Every pair sees the same fixed
oracle, start, acceptance rule, and point-evaluation budget ``40*n+100``.

Large Kostlan cases, including 20x20, use a fixed 3,072-feature surrogate.
They are deliberately not described as exact dense KS systems: exact KS(20,20)
would require C(40,20) = 137,846,528,820 monomials per equation.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import secrets
import statistics
import sys
import time
import warnings
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "benchmarks" / "benchmark_362_blind_adversarial_scipy_root.py"
ENGINE_PATH = ROOT / "flow" / "371_pandrosion_monomial_chain_hirp_standalone.py"
DEFAULT_OUTDIR = ROOT / "benchmarks" / "371_nondeterministic_vs_scipy_root"
METHODS = ["371", "scipy-hybr"]
ACCEPT = 1e-8
MAX_SOLUTION_NORM = 1e8
DEFAULT_PER_FAMILY = 140


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


B = load_module("benchmark371base", BASE_PATH)
E = load_module("pandrosion371benchmark", ENGINE_PATH)

# Reuse the mature checkpointing, SciPy wrapper, CSV writer, and aggregators
# from the sealed 362 harness with the compatible 371 engine.
B.E = E
B.ENGINE_PATH = ENGINE_PATH
B.DEFAULT_OUTDIR = DEFAULT_OUTDIR
B.METHODS = METHODS
B.ACCEPT = ACCEPT
B.MAX_SOLUTION_NORM = MAX_SOLUTION_NORM

BASE_FAMILIES = list(B.FAMILIES)
BASE_MAKE_ORACLE = B.make_oracle
BASE_START_FOR = B.start_for
BASE_ORACLE_FINGERPRINT = B.oracle_fingerprint

EXTRA_FAMILIES = [
    B.Family("monomial_chain_10x4", "monomial-chain", 10, 4, True, "triangular"),
    B.Family("monomial_chain_20x8", "monomial-chain", 20, 8, True, "triangular-flat"),
    B.Family("scaled_chain_12x5", "scaled-chain", 12, 5, True, "equation-spread-1e16"),
    B.Family("permuted_chain_12x5", "permuted-chain", 12, 5, True, "permuted-sparsity"),
    B.Family("ks_feature_20x20_f3072", "feature-kostlan", 20, 20, True, "fixed-bank-surrogate"),
]
FAMILIES = BASE_FAMILIES + EXTRA_FAMILIES
FAMILY_BY_NAME = {family.name: family for family in FAMILIES}
B.FAMILIES = FAMILIES
B.FAMILY_BY_NAME = FAMILY_BY_NAME


class ChainOracle(E.Oracle):
    """Random planted monomial chain, optionally scaled or permuted."""

    def __init__(self, family: Any, seed: int) -> None:
        super().__init__(family.n, family.d, seed, f"blind-{family.kind}-{family.name}")
        self.family = family
        rng = np.random.default_rng(E.mix64(seed + 0x371C4A17))
        phase = rng.uniform(-math.pi, math.pi, family.n)
        self.coefficients = rng.uniform(.7, 1.3, family.n) * np.exp(1j * phase)
        coupling_phase = rng.uniform(-math.pi, math.pi, max(0, family.n - 1))
        self.coupling = rng.uniform(.08, .35, max(0, family.n - 1)) * np.exp(1j * coupling_phase)
        self.upper = bool(rng.integers(0, 2))
        self.permutation = (
            rng.permutation(family.n)
            if family.kind == "permuted-chain"
            else np.arange(family.n)
        )
        self.inverse_permutation = np.argsort(self.permutation)
        internal_root = np.zeros(family.n, np.complex128)
        terminal = family.n - 1 if self.upper else 0
        internal_root[terminal] = rng.uniform(.7, 1.3) * np.exp(
            1j * rng.uniform(-math.pi, math.pi)
        )

        def choose_branch(value: complex) -> complex:
            angle = math.atan2(value.imag, value.real)
            branch = int(rng.integers(0, family.d))
            return abs(value) ** (1.0 / family.d) * np.exp(
                1j * (angle + 2.0 * math.pi * branch) / family.d
            )

        if self.upper:
            for j in range(family.n - 2, -1, -1):
                internal_root[j] = choose_branch(
                    self.coupling[j] * internal_root[j + 1] / self.coefficients[j]
                )
        else:
            for j in range(1, family.n):
                internal_root[j] = choose_branch(
                    self.coupling[j - 1] * internal_root[j - 1] / self.coefficients[j]
                )
        self.internal_root = internal_root
        self.root = np.asarray(internal_root[self.inverse_permutation], np.complex128)
        self.equation_scales = (
            10.0 ** rng.uniform(-8.0, 8.0, family.n)
            if family.kind == "scaled-chain"
            else np.ones(family.n)
        )

    def _unscaled_batch(self, Z: np.ndarray) -> np.ndarray:
        W = np.asarray(Z[:, self.permutation], np.complex128)
        out = self.coefficients[None, :] * W ** self.d
        if self.n > 1:
            if self.upper:
                out[:, :-1] -= self.coupling[None, :] * W[:, 1:]
                out[:, -1] -= self.coefficients[-1] * self.internal_root[-1] ** self.d
            else:
                out[:, 1:] -= self.coupling[None, :] * W[:, :-1]
                out[:, 0] -= self.coefficients[0] * self.internal_root[0] ** self.d
        return np.asarray(out[:, self.inverse_permutation], np.complex128)

    def _eval_batch(self, Z: np.ndarray) -> np.ndarray:
        return self._unscaled_batch(Z) * self.equation_scales[None, :]

    def backward_error(self, z: Any) -> float:
        Z = np.asarray(z, np.complex128)[None, :]
        return E.norm(self._unscaled_batch(Z)[0]) / math.sqrt(max(1, self.n))

    def fingerprint(self) -> str:
        payload = b"".join([
            self.family.name.encode(),
            self.root.view(np.float64).tobytes(),
            self.coefficients.view(np.float64).tobytes(),
            self.coupling.view(np.float64).tobytes(),
            self.permutation.tobytes(),
            self.equation_scales.tobytes(),
            bytes([int(self.upper)]),
        ])
        return B.sha256_bytes(payload)


def make_oracle(problem: Any) -> Any:
    family = FAMILY_BY_NAME[problem.family]
    if family.kind == "feature-kostlan":
        return E.KostlanOracle(
            family.n, family.d, problem.seed, 0, 3072, False, 128
        )
    if family.kind in {"monomial-chain", "scaled-chain", "permuted-chain"}:
        return ChainOracle(family, problem.seed)
    return BASE_MAKE_ORACLE(problem)


def oracle_fingerprint(oracle: Any, family: Any) -> str:
    if isinstance(oracle, ChainOracle):
        return oracle.fingerprint()
    if family.kind == "feature-kostlan":
        payload = (
            family.name.encode()
            + oracle.exps.tobytes()
            + oracle.coeff.view(np.float64).tobytes()
        )
        return B.sha256_bytes(payload)
    return BASE_ORACLE_FINGERPRINT(oracle, family)


def start_for(family: Any, seed: int, index: int) -> np.ndarray:
    if family.kind == "feature-kostlan":
        if index % 4 == 0:
            return np.zeros(family.n, np.complex128)
        radius = (.15, .6, 1.5)[(index - 1) % 3]
        return radius * E.direction(family.n, index, seed)
    if family.kind in {"monomial-chain", "scaled-chain", "permuted-chain"}:
        oracle = ChainOracle(family, seed)
        if index % 3 == 0:
            return np.zeros(family.n, np.complex128)
        radius = (.15, .65)[index % 2]
        return oracle.root + radius * E.direction(family.n, index, seed) / math.sqrt(family.n)
    return BASE_START_FOR(family, seed, index)


B.make_oracle = make_oracle
B.oracle_fingerprint = oracle_fingerprint
B.start_for = start_for


def system_seed(master_seed: int, family_index: int, index: int) -> int:
    mixed = E.mix64(master_seed + 1_000_003 * family_index + 9_176 * index)
    return int(mixed & 0x7FFFFFFF)


def build_manifest(
    per_family: int,
    timing_per_family: int,
    timing_repeats: int,
    master_seed: int,
    seed_source: str,
) -> dict[str, Any]:
    problems: list[dict[str, Any]] = []
    for family_index, family in enumerate(FAMILIES):
        for index in range(per_family):
            seed = system_seed(master_seed, family_index, index)
            start = start_for(family, seed, index)
            shell = B.Problem(
                uid=f"{family.name}_blind{index:04d}",
                family=family.name,
                kind=family.kind,
                n=family.n,
                d=family.d,
                polynomial=family.polynomial,
                condition_class=family.condition_class,
                seed=seed,
                index=index,
                start=tuple(complex(value) for value in start),
                fingerprint="",
            )
            item = asdict(shell)
            item["start"] = B.complex_pairs(shell.start)
            item["fingerprint"] = oracle_fingerprint(make_oracle(shell), family)
            problems.append(item)

    # Randomize campaign order without losing deterministic replay from the
    # recorded entropy seed.  Method order is then alternated by problem.
    rng = np.random.default_rng(E.mix64(master_seed + 0xB11D371))
    rng.shuffle(problems)
    core = {
        "schema": "pandrosion-371-nondeterministic-blind-v1",
        "master_seed": int(master_seed),
        "seed_source": seed_source,
        "per_family": int(per_family),
        "problem_count": len(problems),
        "methods_preregistered": METHODS,
        "hard_budget": "40*n+100 complex point evaluations",
        "acceptance": "finite norm < 1e8 and independent scale-aware backward error < 1e-8",
        "timing_protocol": {
            "repeats": int(timing_repeats),
            "problems_per_family": int(timing_per_family),
            "cluster_unit": "problem",
        },
        "ks_20x20_disclosure": {
            "exact_monomials_per_equation": math.comb(40, 20),
            "backend": "fixed stratified feature bank",
            "features": 3072,
            "exact_dense_ks": False,
            "reason": "exact dense enumeration is computationally infeasible",
        },
        "families": [asdict(family) for family in FAMILIES],
        "problems": problems,
    }
    manifest_hash = B.sha256_bytes(B.canonical_bytes(core))
    return {
        "sealed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "sealed_before_solver_runs": True,
        "manifest_sha256": manifest_hash,
        **core,
    }


def validate_candidate(oracle: Any, z: np.ndarray) -> tuple[float, float, float, bool]:
    solution_norm = E.norm(z)
    raw = E.norm(oracle.eval(z))
    backward = float(oracle.backward_error(z))
    accepted = bool(
        math.isfinite(solution_norm)
        and solution_norm < MAX_SOLUTION_NORM
        and math.isfinite(raw)
        and math.isfinite(backward)
        and backward < ACCEPT
    )
    return solution_norm, raw, backward, accepted


B.validate_candidate = validate_candidate


def engine_args() -> Any:
    args = E.parser().parse_args([])
    args.accept = ACCEPT
    args.validation_accept = ACCEPT
    args.epochs = 96
    args.diagnostic_svd = False
    args.swarm = False
    args.pool = 1
    args.count = 1
    return args


def run_371(problem: Any) -> dict[str, Any]:
    budget = 40 * problem.n + 100
    oracle = make_oracle(problem)
    target = E.Target(oracle)
    args = engine_args()
    start = np.asarray(problem.start, np.complex128)
    caught: list[Any] = []
    started = time.perf_counter()
    exception = None
    budget_exhausted = False
    diagnostics: dict[str, Any] = {}
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("once")
            loc = E.correct(target, start, args, [], args.epochs, None)
        elapsed = time.perf_counter() - started
        solve_evals = int(loc.budget_used)
        solution_norm, raw, backward, accepted = validate_candidate(oracle, loc.y)
        status = loc.status
        library_success = bool(loc.ok)
        iterations = int(loc.epochs)
        telescope = float(loc.max_telescope_defect)
        budget_exhausted = bool(loc.budget_exhausted)
        for name in [
            "halley_proposals", "halley_accepted", "irp_proposals", "irp_accepted",
            "irp2_proposals", "irp2_accepted", "hirp_proposals", "hirp_accepted",
            "monomial_chain_proposals", "monomial_chain_accepted",
            "micro_steps_accepted", "equation_chart_active",
        ]:
            diagnostics[name] = getattr(loc, name, None)
    except Exception as exc:
        elapsed = time.perf_counter() - started
        solve_evals = min(int(oracle.eval_count), budget)
        solution_norm = raw = backward = float("inf")
        accepted = library_success = False
        status, iterations, telescope = "exception", None, None
        exception = repr(exc)
    row = B.common_row(problem, "371", budget)
    row.update({
        "status": status,
        "library_success": library_success,
        "accepted": accepted,
        "seconds": elapsed,
        "solve_evals": solve_evals,
        "solution_norm": solution_norm if math.isfinite(solution_norm) else None,
        "raw_residual": raw if math.isfinite(raw) else None,
        "backward_error": backward if math.isfinite(backward) else None,
        "iterations": iterations,
        "max_telescope_defect": telescope,
        "budget_exceeded": budget_exhausted or solve_evals > budget,
        "warning_count": len(caught),
        "warning_types": ";".join(sorted({w.category.__name__ for w in caught})),
        "exception": exception,
        **diagnostics,
    })
    return row


def run_method(problem: Any, method: str) -> dict[str, Any]:
    if method == "371":
        return run_371(problem)
    return B.run_scipy(problem, method.removeprefix("scipy-"))


B.run_method = run_method


def strict_success(row: dict[str, Any]) -> bool:
    return bool(row["accepted"] and not row["budget_exceeded"])


def paired(
    rows: list[dict[str, Any]], method: str, *, require_clean_termination: bool
) -> dict[str, Any]:
    by_key = {(row["uid"], row["method"]): row for row in rows}
    ids = sorted({row["uid"] for row in rows})
    pairs = [(by_key[(uid, "371")], by_key[(uid, method)]) for uid in ids]
    success = strict_success if require_clean_termination else lambda row: bool(row["accepted"])
    left_only = sum(success(a) and not success(b) for a, b in pairs)
    right_only = sum(success(b) and not success(a) for a, b in pairs)
    discordant = left_only + right_only
    return {
        "paired": len(pairs),
        "both_success": sum(success(a) and success(b) for a, b in pairs),
        "only_371_success": left_only,
        "only_method_success": right_only,
        "neither_success": sum(not success(a) and not success(b) for a, b in pairs),
        "mcnemar_exact_p": float(B.stats.binomtest(left_only, discordant, .5).pvalue) if discordant else 1.0,
    }


def stratified_rate_difference_ci(
    rows: list[dict[str, Any]], *, require_clean_termination: bool, draws: int = 5000
) -> dict[str, Any]:
    """Paired bootstrap, preserving the equal-size family strata."""
    success = strict_success if require_clean_termination else lambda row: bool(row["accepted"])
    by_key = {(row["uid"], row["method"]): row for row in rows}
    groups: list[np.ndarray] = []
    for family in FAMILIES:
        ids = sorted({row["uid"] for row in rows if row["family"] == family.name})
        groups.append(np.asarray([
            float(success(by_key[(uid, "371")]))
            - float(success(by_key[(uid, "scipy-hybr")]))
            for uid in ids
        ]))
    observed = float(np.mean(np.concatenate(groups)))
    rng = np.random.default_rng(0x371B0057 + int(require_clean_termination))
    estimates = np.empty(draws)
    for i in range(draws):
        estimates[i] = float(np.mean(np.concatenate([
            group[rng.integers(0, len(group), len(group))] for group in groups
        ])))
    return {
        "difference": observed,
        "ci95": [float(np.quantile(estimates, .025)), float(np.quantile(estimates, .975))],
        "bootstrap_draws": draws,
        "unit": "paired problem within family strata",
    }


def functional_summary(rows: list[dict[str, Any]], problems: list[Any]) -> dict[str, Any]:
    overall = {
        method: B.aggregate([row for row in rows if row["method"] == method])
        for method in METHODS
    }
    families = []
    for family in FAMILIES:
        subset = [row for row in rows if row["family"] == family.name]
        families.append({
            **asdict(family),
            "problem_count": sum(problem.family == family.name for problem in problems),
            "methods": {
                method: B.aggregate([row for row in subset if row["method"] == method])
                for method in METHODS
            },
        })
    predicates = {
        "polynomial": lambda p: p.polynomial,
        "nonpolynomial": lambda p: not p.polynomial,
        "adversarial": lambda p: p.kind == "adversarial",
        "dense_exact_ks": lambda p: p.kind == "dense-kostlan",
        "fixed_feature_ks": lambda p: p.kind == "feature-kostlan",
        "monomial_chains": lambda p: "chain" in p.kind,
    }
    domains = {}
    for name, predicate in predicates.items():
        uids = {problem.uid for problem in problems if predicate(problem)}
        domains[name] = {
            method: B.aggregate([
                row for row in rows if row["uid"] in uids and row["method"] == method
            ])
            for method in METHODS
        }
    engine_rows = [row for row in rows if row["method"] == "371"]
    counters = [
        "halley_proposals", "halley_accepted", "irp_proposals", "irp_accepted",
        "irp2_proposals", "irp2_accepted", "hirp_proposals", "hirp_accepted",
        "monomial_chain_proposals", "monomial_chain_accepted", "micro_steps_accepted",
    ]
    diagnostics = {
        name: sum(int(row.get(name) or 0) for row in engine_rows) for name in counters
    }
    diagnostics["equation_chart_active_runs"] = sum(
        bool(row.get("equation_chart_active")) for row in engine_rows
    )
    return {
        "overall": overall,
        "families": families,
        "domains": domains,
        "paired_vs_371": {
            method: paired(rows, method, require_clean_termination=True)
            for method in METHODS[1:]
        },
        "paired_validated_vs_371": {
            method: paired(rows, method, require_clean_termination=False)
            for method in METHODS[1:]
        },
        "rate_difference_ci95": {
            "validated": stratified_rate_difference_ci(
                rows, require_clean_termination=False
            ),
            "strict": stratified_rate_difference_ci(
                rows, require_clean_termination=True
            ),
        },
        "371_diagnostics": diagnostics,
    }


def timing_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    per_problem: dict[tuple[str, str], list[float]] = {}
    per_problem_success: dict[tuple[str, str], list[bool]] = {}
    for row in rows:
        key = (row["uid"], row["method"])
        per_problem.setdefault(key, []).append(float(row["seconds"]))
        per_problem_success.setdefault(key, []).append(strict_success(row))
    medians = {key: statistics.median(values) for key, values in per_problem.items()}
    stable = {key: all(values) for key, values in per_problem_success.items()}
    rng = np.random.default_rng(0x371C1)
    ids = sorted({uid for uid, _ in medians})
    methods = {}
    for method in METHODS:
        values = [medians[(uid, method)] for uid in ids if (uid, method) in medians]
        successful = [
            medians[(uid, method)] for uid in ids if stable.get((uid, method), False)
        ]
        common_ids = [] if method == "371" else [
            uid for uid in ids
            if stable.get((uid, "371"), False)
            and stable.get((uid, method), False)
            and medians[(uid, method)] > 0
        ]
        ratios = [medians[(uid, "371")] / medians[(uid, method)] for uid in common_ids]
        methods[method] = {
            "problem_count": len(values),
            "repeat_count": len([row for row in rows if row["method"] == method]),
            "stable_success_problem_count": len(successful),
            "common_success_problem_count": None if method == "371" else len(common_ids),
            "median_termination_seconds_ci95": B.bootstrap_median_ci(values, rng),
            "median_stable_success_seconds_ci95": B.bootstrap_median_ci(successful, rng) if successful else None,
            "paired_ratio_371_over_method_common_success_ci95": B.bootstrap_median_ci(ratios, rng) if ratios else None,
        }
    return {"cluster_unit": "problem", "bootstrap_draws": 5000, "methods": methods}


def fmt(value: Any) -> str:
    if value is None:
        return "—"
    x = float(value)
    return f"{x:.2e}" if x and (abs(x) < 1e-3 or abs(x) >= 1e4) else f"{x:.3f}"


def render_report(payload: dict[str, Any]) -> str:
    summary = payload["functional_summary"]
    timing = payload["timing_summary"]
    meta = payload["metadata"]
    manifest = payload["manifest"]
    lines = [
        "# Pandrosion 371 vs `scipy.optimize.root(hybr)` — benchmark non déterministe", "",
        f"Date : {meta['timestamp']}", "",
        "## Protocole", "",
        f"Le manifeste de **{meta['problem_count']} systèmes** a été généré depuis l'entropie système, "
        "écrit et haché avant la première exécution.",
        f"Graine maître rejouable : `{manifest['master_seed']}`.",
        f"SHA-256 du manifeste : `{meta['manifest_sha256']}`.", "",
        f"{len(manifest['families'])} familles × {manifest['per_family']} réalisations. Même système, même départ, "
        "même plafond `40*n+100` évaluations de points complexes et même validation externe pour les deux solveurs.",
        "SciPy reçoit la réalification C^n → R^(2n), sans Jacobien fourni. Le succès publié dépend de l'erreur "
        "arrière indépendante, pas du drapeau interne du solveur.", "",
        "KS 20×20 est une **banque fixe de 3 072 caractéristiques**, identique pour les deux solveurs. "
        f"Ce n'est pas un KS dense exact : celui-ci demanderait {manifest['ks_20x20_disclosure']['exact_monomials_per_equation']:,} "
        "monômes par équation.", "",
        "## Résultats globaux", "",
        "| méthode | racines validées | succès stricts | taux strict | évaluations médianes | p95 évaluations | temps médian | plafonds | exceptions |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for method in METHODS:
        item = summary["overall"][method]
        lines.append(
            f"| {method} | {item['accepted']}/{item['attempted']} | {item['strict_successes']}/{item['attempted']} "
            f"| {100 * item['strict_success_rate']:.1f} % "
            f"| {fmt(item['median_evals'])} | {fmt(item['p95_evals'])} | {fmt(item['median_seconds'])} s "
            f"| {item['budget_hits']} | {item['exceptions']} |"
        )
    validated_effect = summary["rate_difference_ci95"]["validated"]
    strict_effect = summary["rate_difference_ci95"]["strict"]
    validated_pair = summary["paired_validated_vs_371"]["scipy-hybr"]
    lines += [
        "",
        "Une racine validée a été observée dans le budget pour "
        f"{summary['overall']['371']['accepted']}/{summary['overall']['371']['attempted']} cas avec 371 et "
        f"{summary['overall']['scipy-hybr']['accepted']}/{summary['overall']['scipy-hybr']['attempted']} avec SciPy. "
        f"Écart apparié : **{100 * validated_effect['difference']:.1f} points** "
        f"(IC95 bootstrap stratifié [{100 * validated_effect['ci95'][0]:.1f}, {100 * validated_effect['ci95'][1]:.1f}]).",
        f"Décomposition validée : {validated_pair['both_success']} succès communs, "
        f"{validated_pair['only_371_success']} pour 371 seul, {validated_pair['only_method_success']} pour SciPy seul, "
        f"{validated_pair['neither_success']} échecs communs. La version conservatrice exigeant une terminaison "
        f"sans signal de plafond donne un écart de {100 * strict_effect['difference']:.1f} points "
        f"(IC95 [{100 * strict_effect['ci95'][0]:.1f}, {100 * strict_effect['ci95'][1]:.1f}]).",
    ]
    lines += [
        "", "## Comparaison appariée", "",
        "Le test exact de McNemar porte sur les succès stricts du même système et du même point initial.", "",
        "| méthode SciPy | succès communs | 371 seul | SciPy seule | échecs communs | p exacte |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for method in METHODS[1:]:
        item = summary["paired_vs_371"][method]
        p_text = "<5e-324" if item["mcnemar_exact_p"] == 0 else fmt(item["mcnemar_exact_p"])
        lines.append(
            f"| {method} | {item['both_success']} | {item['only_371_success']} | "
            f"{item['only_method_success']} | {item['neither_success']} | {p_text} |"
        )
    lines += [
        "", "## Résultats par domaine", "",
        "| domaine | méthode | succès stricts | taux |",
        "|---|---|---:|---:|",
    ]
    for domain, methods in summary["domains"].items():
        for method in METHODS:
            item = methods[method]
            rate = 100 * item["strict_success_rate"] if item["attempted"] else 0.0
            lines.append(
                f"| {domain} | {method} | {item['strict_successes']}/{item['attempted']} | {rate:.1f} % |"
            )
    lines += [
        "", "## Détail par famille", "",
        "| famille | n | degré | nature | 371 | SciPy hybr |",
        "|---|---:|---:|---|---:|---:|",
    ]
    for family in summary["families"]:
        total = family["problem_count"]
        a = family["methods"]["371"]["strict_successes"]
        b = family["methods"]["scipy-hybr"]["strict_successes"]
        lines.append(
            f"| {family['name']} | {family['n']} | {family['d']} | {family['condition_class']} | {a}/{total} | {b}/{total} |"
        )
    ks_family = next(
        family for family in summary["families"]
        if family["name"] == "ks_feature_20x20_f3072"
    )
    lines += [
        "", "## Focus KS-feature 20×20", "",
        "| méthode | succès stricts | évaluations médianes | p95 évaluations | temps répété médian |",
        "|---|---:|---:|---:|---:|",
    ]
    timing_runs = [
        row for row in payload.get("timing_runs", [])
        if row["family"] == "ks_feature_20x20_f3072"
    ]
    for method in METHODS:
        item = ks_family["methods"][method]
        by_uid: dict[str, list[float]] = {}
        for row in timing_runs:
            if row["method"] == method:
                by_uid.setdefault(row["uid"], []).append(float(row["seconds"]))
        repeated = [statistics.median(values) for values in by_uid.values()]
        repeated_median = statistics.median(repeated) if repeated else None
        lines.append(
            f"| {method} | {item['strict_successes']}/{item['attempted']} | "
            f"{fmt(item['median_evals'])} | {fmt(item['p95_evals'])} | {fmt(repeated_median)} s |"
        )
    if timing.get("methods"):
        lines += [
            "", "## Chronométrage répété", "",
            "IC 95 % par bootstrap de 5 000 tirages sur les médianes par problème.", "",
            "| méthode | terminaison médiane [IC95] s | succès stables | ratio 371/méthode [IC95] | n commun |",
            "|---|---:|---:|---:|---:|",
        ]
        for method in METHODS:
            item = timing["methods"][method]
            ci = item["median_termination_seconds_ci95"]
            ratio = item["paired_ratio_371_over_method_common_success_ci95"]
            ratio_text = "—" if ratio is None else f"{fmt(ratio[0])} [{fmt(ratio[1])}, {fmt(ratio[2])}]"
            common = "—" if item["common_success_problem_count"] is None else str(item["common_success_problem_count"])
            lines.append(
                f"| {method} | {fmt(ci[0])} [{fmt(ci[1])}, {fmt(ci[2])}] | "
                f"{item['stable_success_problem_count']}/{item['problem_count']} | {ratio_text} | {common} |"
            )
    diag = summary["371_diagnostics"]
    lines += [
        "", "## Diagnostics 371", "",
        f"Charts d'équations activés : {diag['equation_chart_active_runs']} exécutions. "
        f"Chaînes monomiales acceptées : {diag['monomial_chain_accepted']}/{diag['monomial_chain_proposals']}. "
        f"HIRP acceptés : {diag['hirp_accepted']}/{diag['hirp_proposals']}. "
        f"IRP-2 acceptés : {diag['irp2_accepted']}/{diag['irp2_proposals']}.", "",
        "## Limites", "",
        "Le corpus est généré par ce harnais, pas par un tiers indépendant. Le caractère non déterministe vient "
        "de la graine d'entropie fraîche ; le manifeste scellé rend néanmoins cette campagne exactement rejouable. "
        "Les grandes cases KS évaluent une approximation à banque fixe et non le polynôme KS dense exact. "
        "Les temps muraux sont propres à cette machine ; les comptes d'évaluations sont la mesure portable principale.", "",
        f"SHA-256 du moteur 371 : `{meta['engine_sha256']}`.", "",
        "Artefacts : `blind_manifest.json`, `functional_raw.csv`, `timing_raw.csv`, `benchmark.json`, `report.md`.", "",
    ]
    return "\n".join(lines)


def validate_results(rows: list[dict[str, Any]], problems: list[Any]) -> None:
    expected = {(problem.uid, method) for problem in problems for method in METHODS}
    observed = [(row["uid"], row["method"]) for row in rows]
    if len(observed) != len(expected) or set(observed) != expected:
        raise AssertionError("incomplete or duplicate functional results")
    for row in rows:
        if row["solve_evals"] > row["eval_budget"]:
            raise AssertionError(f"hard budget crossed by {row['uid']} {row['method']}")
        if row["accepted"]:
            if row["solution_norm"] is None or row["solution_norm"] >= MAX_SOLUTION_NORM:
                raise AssertionError("invalid accepted solution norm")
            if row["backward_error"] is None or row["backward_error"] >= ACCEPT:
                raise AssertionError("invalid accepted backward error")


def write_overall_csv(path: Path, summary: dict[str, Any]) -> None:
    rows = [{"method": method, **summary["overall"][method]} for method in METHODS]
    B.write_csv(path, rows)


def parse_seed(raw: str) -> int:
    value = int(raw, 0)
    if value < 0 or value >= 2**64:
        raise argparse.ArgumentTypeError("seed must be in [0, 2**64)")
    return value


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--per-family", type=int, default=DEFAULT_PER_FAMILY)
    parser.add_argument("--create-manifest", action="store_true")
    parser.add_argument("--seed", type=parse_seed)
    parser.add_argument("--functional-only", action="store_true")
    parser.add_argument("--timing-only", action="store_true")
    parser.add_argument("--timing-per-family", type=int, default=10)
    parser.add_argument("--timing-repeats", type=int, default=7)
    parser.add_argument("--progress-every", type=int, default=25)
    args = parser.parse_args(argv)
    if args.per_family <= 0 or args.timing_per_family <= 0 or args.timing_repeats <= 0:
        parser.error("all benchmark sizes must be positive")
    if args.functional_only and args.timing_only:
        parser.error("--functional-only and --timing-only are mutually exclusive")
    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.outdir / "blind_manifest.json"

    if args.create_manifest:
        if manifest_path.exists():
            raise FileExistsError(f"refusing to overwrite sealed manifest {manifest_path}")
        master_seed = args.seed if args.seed is not None else secrets.randbits(64)
        seed_source = "user" if args.seed is not None else "os-entropy"
        manifest = build_manifest(
            args.per_family, args.timing_per_family, args.timing_repeats,
            master_seed, seed_source,
        )
        manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False), encoding="utf-8")
        print(
            f"sealed {manifest['problem_count']} problems seed={master_seed} "
            f"sha256={manifest['manifest_sha256']}",
            flush=True,
        )
        return 0

    manifest, problems = B.load_manifest(manifest_path)
    if manifest["methods_preregistered"] != METHODS:
        raise ValueError("methods differ from the sealed protocol")
    if [item["name"] for item in manifest["families"]] != [item.name for item in FAMILIES]:
        raise ValueError("families differ from the sealed protocol")
    protocol = manifest["timing_protocol"]
    if (
        args.timing_per_family != protocol["problems_per_family"]
        or args.timing_repeats != protocol["repeats"]
    ):
        raise ValueError("timing arguments differ from the sealed protocol")

    started = time.perf_counter()
    functional_rows = B.read_jsonl(args.outdir / "functional_results.jsonl")
    if not args.timing_only:
        functional_rows = B.run_functional(problems, args.outdir, args.progress_every)
        validate_results(functional_rows, problems)
        B.write_csv(args.outdir / "functional_raw.csv", functional_rows)
    elif not functional_rows:
        raise FileNotFoundError("functional checkpoint required for --timing-only")

    timing_problems = B.timing_selection(problems, args.timing_per_family)
    timing_rows = B.read_jsonl(args.outdir / "timing_results.jsonl")
    if not args.functional_only:
        timing_rows = B.run_timing(
            timing_problems, args.outdir, args.timing_repeats, args.progress_every
        )
        B.validate_timing(timing_rows, timing_problems, args.timing_repeats)
        B.write_csv(args.outdir / "timing_raw.csv", timing_rows)

    functional = functional_summary(functional_rows, problems)
    timing = timing_summary(timing_rows) if timing_rows else {"methods": {}}
    metadata = {
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "platform": platform.platform(),
        "problem_count": len(problems),
        "solver_run_count": len(functional_rows),
        "timing_problem_count": len(timing_problems),
        "timing_run_count": len(timing_rows),
        "manifest_sha256": manifest["manifest_sha256"],
        "master_seed": manifest["master_seed"],
        "engine_sha256": hashlib.sha256(ENGINE_PATH.read_bytes()).hexdigest(),
        "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "benchmark_seconds_this_invocation": time.perf_counter() - started,
    }
    public_manifest = {key: value for key, value in manifest.items() if key != "problems"}
    payload = {
        "metadata": metadata,
        "manifest": public_manifest,
        "functional_summary": functional,
        "timing_summary": timing,
        "functional_runs": functional_rows,
        "timing_runs": timing_rows,
    }
    (args.outdir / "benchmark.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8"
    )
    (args.outdir / "report.md").write_text(render_report(payload), encoding="utf-8")
    write_overall_csv(args.outdir / "summary.csv", functional)
    print(f"wrote {args.outdir / 'report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
