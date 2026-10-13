#!/usr/bin/env python3
"""Reproducible DD convergence diagnostics with LinearAlifold enabled.

No production edits. Uses the existing runner and paired-summary helpers.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from statistics import fmean

from benchmark_dd_budget import PILOT_DATASET_IDS, paired_summary as original_paired_summary, write_csv
from benchmark_luna_optimization import trace_metrics

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "benchmarks/configs/murlet-dataset1-linearalifold-repair-v2-iter100-beam100.json"
RESULTS = ROOT / "benchmarks/results-linearalifold-dd-convergence-20260919"
BASES = (
    "cl_linear_dd_0_alifold", "lpc_contralign_w2_alifold",
    "lpv_linear_dd_0p01_alifold", "lpv_probcons_w2_alifold",
)
BUDGETS = (100, 200, 300, 500)
CHECKPOINTS = (1, 10, 25, 50, 100, 200, 300, 500)
ITERATION_HEADER = re.compile(r'"merge_id":(\d+),"iteration":(\d+),')


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def paired_summary(candidate, reference, label):
    # The older helper's per-condition subsection has its former experiment's
    # four IDs hardcoded. Each comparison here already covers exactly one base.
    assert len({r["condition"] for r in candidate}) == 1
    result = original_paired_summary(candidate, reference, label)
    result.pop("per_condition", None)
    assert result["paired_count"] == len(candidate) == len(reference)
    assert result["missing_candidate"] == result["missing_reference"] == 0
    return result


def prepare():
    if RESULTS.exists():
        raise RuntimeError(f"Refusing to overwrite {RESULTS}")
    RESULTS.mkdir(parents=True)
    binary = RESULTS / "dafs-frozen"
    shutil.copy2(ROOT / "build-arm64/src/dafs", binary)
    source = read(SOURCE)
    original = {c["id"]: c for c in source["conditions"]}
    inputs = ROOT / "benchmarks/data/murlet/manifest.json"
    assert sum(d.get("collection") == "dataset1" for d in read(inputs)["datasets"]) == 85
    common = {
        "binary": str(binary), "timeout_seconds": 3600,
        "repetitions": 1, "random_seed": 42,
        "environment": {"OMP_NUM_THREADS": "1"},
        "sci": source["sci"], "dataset_manifests": [str(inputs)],
    }
    for scope in ("alifold", "matched-off-pilot"):
        config = copy.deepcopy(common)
        config["output_dir"] = str(RESULTS / scope)
        config["dataset_filter"] = ({"collections": ["dataset1"]}
                                    if scope == "alifold" else {"ids": list(PILOT_DATASET_IDS)})
        config["conditions"] = []
        for base in BASES:
            for mode in (("normal", "projected") if scope == "alifold" else ("normal",)):
                for budget in (BUDGETS if scope == "alifold" else (100, 500)):
                    c = copy.deepcopy(original[base])
                    args = [f"--max-iter={budget}" if a.startswith("--max-iter=") else a
                            for a in c["args"]]
                    if scope != "alifold":
                        args = ["--no-alifold" if a == "--alifold" else a for a in args]
                    if mode == "projected":
                        args.append("--dd-projected-norm")
                    args.extend(["--seed=42", "--refinement=0"])
                    c.update(id=f"{base}_{mode}_{budget}", args=args,
                             base_condition=base, mode=mode, max_iter=budget)
                    config["conditions"].append(c)
        config["convergence_experiment"] = {
            "source_config": str(SOURCE), "source_sha256": digest(SOURCE),
            "binary_sha256": digest(binary), "scope": scope,
            "expected_runs": len(config["conditions"]) * (85 if scope == "alifold" else 12),
            "production_changes": False, "beam": 100,
            "fixed_stepsize": 0.5, "bp_update": False,
        }
        save(RESULTS / f"{scope}-config.json", config)
    paths = [p for p in (ROOT / "src").rglob("*") if p.suffix in (".cpp", ".h", ".hpp")]
    paths += [ROOT / "CMakeLists.txt", ROOT / "benchmarks/run.py", SOURCE, inputs]
    paths += list((ROOT / "src").rglob("CMakeLists.txt"))
    save(RESULTS / "source-hashes.json", {str(p.relative_to(ROOT)): digest(p) for p in sorted(set(paths))})
    print(json.dumps({"binary_sha256": digest(binary), "full_runs": 2720,
                      "matched_off_control_runs": 96, "root": str(RESULTS)}))


def rows_for(scope):
    config = read(RESULTS / f"{scope}-config.json")
    manifest = read(RESULTS / scope / "manifest.json")
    assert manifest["binary_sha256"] == config["convergence_experiment"]["binary_sha256"] == digest(RESULTS / "dafs-frozen")
    conditions = {c["id"]: c for c in config["conditions"]}
    rows, checkpoints = [], []
    for path in sorted((RESULTS / scope / "runs").glob("*/*/rep-*/result.json")):
        result = read(path)
        c = conditions[result["condition"]]
        summaries, last_lines, run_summary, stages = [], {}, {}, {}
        def checkpoint(e, iteration, carried=False):
            return {
                "scope": scope, "condition": c["base_condition"], "mode": c["mode"],
                "dataset": result["dataset"], "merge": e["merge_id"],
                "iteration": iteration, "observed_iteration": e["iteration"] + 1,
                "carried_after_stop": carried, "lb": e["lb"], "ub": e["best_ub"],
                "gap": e["gap"], "relative_gap": e["gap"] / max(1, abs(e["lb"])),
                "violations": e["violated"], "certified_violations": e["certified_track_violated"],
            }
        with (path.parent / "metrics.jsonl").open() as stream:
            for line in stream:
                # Runner already validates every event. Only materialize the
                # terminal/checkpoint events needed here, avoiding millions
                # of otherwise discarded JSON dictionaries.
                header = ITERATION_HEADER.search(line) if '"event":"dd_iteration"' in line else None
                if header is not None:
                    merge_id, iteration = map(int, header.groups())
                    last_lines[merge_id] = line
                    if c["max_iter"] == 500 and iteration + 1 in CHECKPOINTS:
                        checkpoints.append(checkpoint(json.loads(line), iteration + 1))
                    continue
                e = json.loads(line)
                event = e["event"]
                if event == "dd_summary":
                    summaries.append(e)
                elif event == "dd_iteration":
                    last_lines[e["merge_id"]] = line
                    if c["max_iter"] == 500 and e["iteration"] + 1 in CHECKPOINTS:
                        checkpoints.append(checkpoint(e, e["iteration"] + 1))
                elif event == "run_summary":
                    run_summary = e
                elif event == "stage":
                    stages[e.get("name")] = e["seconds"]
        last = {key: json.loads(line) for key, line in last_lines.items()}
        # Fixed cohort: early-stopped merges retain their terminal observation.
        # This avoids comparing progressively harder surviving merges at later t.
        if c["max_iter"] == 500:
            for e in last.values():
                for iteration in CHECKPOINTS:
                    if iteration > e["iteration"] + 1:
                        checkpoints.append(checkpoint(e, iteration, carried=True))
        metrics = trace_metrics([*summaries, *last.values(), run_summary], result)
        rel = [e["gap"] / max(1, abs(e["lb"])) for e in summaries]
        reasons = Counter(e["stop_reason"] for e in summaries)
        row = {
            "scope": scope, "condition": c["base_condition"], "mode": c["mode"],
            "budget": c["max_iter"], "dataset": result["dataset"], "repetition": result["repetition"],
            "status": result["status"], "metrics_valid": result["metrics_validation"]["valid"],
            "quality_complete": result["quality_complete"], **metrics,
            "relative_gap_mean": fmean(rel), "relative_gap_max": max(rel),
            "certified_stop_merges": reasons["certified_gap"],
            "max_iter_merges": reasons["max_iterations"], "stationary_merges": reasons["beam_stationary"],
            "initial_cbp_mean": fmean(e["initial_cbp"] for e in summaries),
            "final_cbp_mean": fmean(e["final_cbp"] for e in summaries),
            "cbp_peak": run_summary.get("cbp_peak"),
            "bpp_seconds": stages.get("base_pair_probabilities"),
            "progressive_seconds": stages.get("progressive_alignment"),
            "final_structure_seconds": stages.get("final_consensus_folding"),
        }
        rows.append(row)
    return rows, checkpoints


def summarize():
    rows, checkpoints = [], []
    for scope in ("alifold", "matched-off-pilot"):
        r, c = rows_for(scope)
        rows.extend(r)
        checkpoints.extend(c)
    assert len(rows) == 2816, len(rows)
    assert len({(r["scope"], r["condition"], r["mode"], r["budget"], r["dataset"]) for r in rows}) == len(rows)
    if not all(r["status"] == "success" and r["metrics_valid"] and r["quality_complete"] for r in rows):
        raise RuntimeError("Invalid/incomplete benchmark results; inspect raw results")
    groups = defaultdict(list)
    for r in rows:
        groups[r["scope"], r["condition"], r["mode"], r["budget"]].append(r)
    assert len(groups) == 40
    assert all(len(v) == (85 if k[0] == "alifold" else 12) for k, v in groups.items())
    metrics = ("mcc", "sci", "sps", "wall_seconds", "cpu_seconds", "dd_wall_seconds", "dd_iterations",
               "dd_violations_total", "certified_violations_total", "dd_final_gap_mean", "dd_final_gap_last",
               "relative_gap_mean", "certified_stop_merges", "max_iter_merges", "stationary_merges",
               "initial_cbp_mean", "final_cbp_mean", "cbp_peak",
               "bpp_seconds", "progressive_seconds", "final_structure_seconds")
    summary = []
    for key, items in sorted(groups.items()):
        record = dict(zip(("scope", "condition", "mode", "budget"), key), n=len(items))
        for m in metrics:
            vals = [r[m] for r in items if r.get(m) is not None]
            record[m] = fmean(vals) if vals else None
        summary.append(record)
    comparisons = []
    for base in BASES:
        for mode in ("normal", "projected"):
            for budget in BUDGETS:
                candidate = groups["alifold", base, mode, budget]
                comparisons.append(paired_summary(candidate, groups["alifold", base, "normal", 500],
                                                   f"{base}:{mode}:{budget}:vs-normal500"))
                if mode == "projected":
                    comparisons.append(paired_summary(candidate, groups["alifold", base, "normal", budget],
                                                       f"{base}:{mode}:{budget}:vs-same-budget"))
        for budget in (100, 500):
            candidate = [r for r in groups["alifold", base, "normal", budget] if r["dataset"] in PILOT_DATASET_IDS]
            comparisons.append(paired_summary(candidate, groups["matched-off-pilot", base, "normal", budget],
                                               f"{base}:normal:{budget}:alifold-vs-off-matched12"))
    changed = [p for p, sha in read(RESULTS / "source-hashes.json").items() if digest(ROOT / p) != sha]
    report = {"runs": len(rows), "valid_runs": len(rows), "summary": summary, "comparisons": comparisons,
              "source_files_changed_since_freeze": changed,
              "binary_sha256": digest(RESULTS / "dafs-frozen"),
              "analysis_hashes": {name: digest(ROOT / "benchmarks" / name) for name in
                                  ("benchmark_linearalifold_dd.py", "benchmark_dd_budget.py", "benchmark_luna_optimization.py")},
              "scope_note": "Full85 alifold sweep; matched off controls only12 cases, not old tuned-off parameters. No convergence theorem; objective fixed within DD merge only."}
    save(RESULTS / "report.json", report)
    write_csv(RESULTS / "summary-convergence.csv", summary)
    write_csv(RESULTS / "case-convergence.csv", rows)
    write_csv(RESULTS / "within-merge-checkpoints.csv", checkpoints)
    curve_groups = defaultdict(list)
    for r in checkpoints:
        curve_groups[r["scope"], r["condition"], r["mode"], r["iteration"]].append(r)
    curves = []
    for key, items in sorted(curve_groups.items()):
        entry = dict(zip(("scope", "condition", "mode", "iteration"), key), merges=len(items))
        for metric in ("gap", "relative_gap", "violations", "certified_violations"):
            entry[metric] = fmean(r[metric] for r in items)
        entry["carried_after_stop"] = sum(r["carried_after_stop"] for r in items)
        curves.append(entry)
    write_csv(RESULTS / "within-merge-curves.csv", curves)
    print(json.dumps({"runs": len(rows), "comparisons": len(comparisons), "changed_sources": changed}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "summarize"))
    args = parser.parse_args()
    prepare() if args.command == "prepare" else summarize()
