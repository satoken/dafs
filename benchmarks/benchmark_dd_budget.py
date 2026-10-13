#!/usr/bin/env python3
"""Generate and summarize the DD max-iteration budget experiment.

Bulk execution is deliberately delegated to the existing benchmarks/run.py
runner. This driver only materializes experiment configs and reads completed
result trees. It does not build or modify production sources.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from statistics import fmean, median
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"
SOURCE_CONFIG = BENCHMARKS / "configs/luna-relaxed.json"
DEFAULT_BINARY = ROOT / "build-arm64/src/dafs"
DEFAULT_RESULTS_ROOT = BENCHMARKS / "results-murlet-dataset1-dd-budget-arm64-v272"
EXPECTED_BINARY_SHA256 = (
    "753129a65753f63b33031c3284ce00be5a09559706c7c21955e6d8d0df2e1695"
)

MODES = ("baseline", "projected-norm", "column-bound", "tight-alignment")
MODE_FLAGS = {
    "baseline": None,
    "projected-norm": "--dd-projected-norm",
    "column-bound": "--dd-column-bound",
    "tight-alignment": "--dd-tight-alignment",
}
MODE_COMPLEXITY = {
    "baseline": "existing default path",
    "projected-norm": "strict-linear candidate: O(L+M) accounting path",
    "column-bound": "strict-linear candidate: O(L+M) column-capacity path",
    "tight-alignment": "reference only: O(M log L) exact sparse monotone alignment bound",
}
BUDGETS = (100, 200, 300, 500)
BASE_CONDITIONS = (
    "cl_linear_ath02",
    "lpv_probcons_fpct05",
    "lpc_contralign_ath02",
    "lpv_linear_fpct05",
)

# Frozen before inspecting this experiment's quality results. These are the
# twelve cases used by the preceding DD pilot, spanning families and lengths.
PILOT_DATASET_IDS = (
    "murlet-dataset1-RF00436-2",
    "murlet-dataset1-RF00031-2",
    "murlet-dataset1-RF00181-2",
    "murlet-dataset1-RF00169-2",
    "murlet-dataset1-RF00059-2",
    "murlet-dataset1-RF00001-2",
    "murlet-dataset1-RF00175-2",
    "murlet-dataset1-RF00050-2",
    "murlet-dataset1-RF00002-2",
    "murlet-dataset1-RF00004-2",
    "murlet-dataset1-RF00230-1",
    "murlet-dataset1-RF00017-1",
)

QUALITY_FIELDS = ("mcc", "sci", "sps")
SUMMARY_FIELDS = (
    *QUALITY_FIELDS,
    "dd_iterations",
    "dd_wall_seconds",
    "wall_seconds",
    "cpu_seconds",
    "max_rss_kb",
    "dd_violations_total",
    "dd_violations_mean_merge",
    "dd_violations_max_merge",
    "dd_final_gap_mean",
    "dd_final_gap_last",
    "dd_final_gap_max",
    "certified_violations_total",
    "certified_violations_mean_merge",
    "certified_final_gap_mean",
    "certified_final_gap_max",
)
FLAG_NAMES = tuple(value for value in MODE_FLAGS.values() if value)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any, force: bool = True) -> None:
    if path.exists() and not force:
        raise RuntimeError(f"refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def absolute(value: str | Path, base: Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def normalized_source(path: Path) -> dict[str, Any]:
    config = copy.deepcopy(read_json(path))
    base = path.parent
    config["dataset_manifests"] = [
        str(absolute(item, base)) for item in config.get("dataset_manifests", [])
    ]
    if config.get("sci", {}).get("bin_dir"):
        config["sci"]["bin_dir"] = str(absolute(config["sci"]["bin_dir"], base))
    return config


def all_datasets(config: dict[str, Any]) -> list[dict[str, Any]]:
    result = copy.deepcopy(config.get("datasets", []))
    for manifest_name in config.get("dataset_manifests", []):
        manifest = read_json(Path(manifest_name))
        result.extend(copy.deepcopy(manifest["datasets"]))
    selected = config.get("dataset_filter", {})
    if selected.get("collections"):
        allowed = {str(item) for item in selected["collections"]}
        result = [item for item in result if str(item.get("collection")) in allowed]
    if selected.get("ids"):
        allowed = {str(item) for item in selected["ids"]}
        result = [item for item in result if str(item.get("id")) in allowed]
    return result


def max_iter(condition: dict[str, Any]) -> int | None:
    if condition.get("max_iter") is not None:
        return int(condition["max_iter"])
    args = [str(item) for item in condition.get("args", [])]
    for index, argument in enumerate(args):
        match = re.fullmatch(r"--max-iter(?:=(\d+))?", argument)
        if match:
            if match.group(1):
                return int(match.group(1))
            if index + 1 < len(args):
                return int(args[index + 1])
    return None


def replace_max_iter(args: list[str], budget: int) -> list[str]:
    result = list(args)
    for index, argument in enumerate(result):
        if argument.startswith("--max-iter="):
            result[index] = f"--max-iter={budget}"
            break
        if argument == "--max-iter":
            if index + 1 >= len(result):
                raise ValueError("--max-iter has no value")
            result[index + 1] = str(budget)
            break
    else:
        raise ValueError("source condition has no --max-iter")
    return result


def remove_option(
    args: list[str], boolean_names: set[str], value_names: set[str]
) -> list[str]:
    result: list[str] = []
    index = 0
    while index < len(args):
        argument = str(args[index])
        if argument in boolean_names:
            index += 1
            continue
        if argument in value_names:
            index += 2
            continue
        if any(argument.startswith(f"{name}=") for name in value_names):
            index += 1
            continue
        result.append(argument)
        index += 1
    return result


def condition_for(
    source_condition: dict[str, Any], mode: str, budget: int
) -> dict[str, Any]:
    condition = copy.deepcopy(source_condition)
    args = replace_max_iter([str(item) for item in condition.get("args", [])], budget)
    args = remove_option(args, set(FLAG_NAMES), {"--seed", "--refinement"})
    args.extend(("--seed=42", "--refinement=0"))
    flag = MODE_FLAGS[mode]
    if flag:
        args.append(flag)
    base = str(condition["base_condition"])
    condition["id"] = f"{base}_iter{budget}"
    condition["args"] = args
    condition["max_iter"] = budget
    condition["dd_budget_mode"] = mode
    return condition


def source_conditions(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result = {}
    for condition in config.get("conditions", []):
        base = str(condition.get("base_condition", condition["id"]))
        if base in BASE_CONDITIONS:
            result[base] = condition
    missing = [item for item in BASE_CONDITIONS if item not in result]
    if missing:
        raise ValueError(f"source config is missing conditions: {missing}")
    return result


def plan(binary: Path) -> dict[str, Any]:
    source = normalized_source(SOURCE_CONFIG)
    datasets = all_datasets(source)
    dataset_ids = [str(item["id"]) for item in datasets]
    missing = [item for item in PILOT_DATASET_IDS if item not in dataset_ids]
    if missing:
        raise ValueError(f"pilot datasets missing from dataset1: {missing}")
    return {
        "schema_version": 1,
        "experiment": "murlet-dataset1-dd-budget",
        "source_config": str(SOURCE_CONFIG.resolve()),
        "source_config_sha256": sha256_file(SOURCE_CONFIG),
        "binary": str(binary.resolve()),
        "binary_sha256": sha256_file(binary) if binary.is_file() else None,
        "expected_binary_sha256": EXPECTED_BINARY_SHA256,
        "dataset_collection": "dataset1",
        "dataset_count": len(dataset_ids),
        "pilot_dataset_count": len(PILOT_DATASET_IDS),
        "base_conditions": list(BASE_CONDITIONS),
        "budgets": list(BUDGETS),
        "modes": {
            mode: {"flag": MODE_FLAGS[mode], "complexity_note": MODE_COMPLEXITY[mode]}
            for mode in MODES
        },
        "pilot_runs": len(PILOT_DATASET_IDS) * len(BASE_CONDITIONS) * len(BUDGETS) * len(MODES),
        "full_baseline_runs": len(dataset_ids) * len(BASE_CONDITIONS) * 2,
        "full_candidate_runs_per_mode_budget": len(dataset_ids) * len(BASE_CONDITIONS),
        "full_followup_runs": len(dataset_ids) * 4,
        "full_followup_arms": [
            {
                "mode": "projected-norm",
                "budget": 200,
                "base_condition": "lpc_contralign_ath02",
            },
            {
                "mode": "projected-norm",
                "budget": 300,
                "base_condition": "lpv_probcons_fpct05",
            },
            {
                "mode": "baseline",
                "budget": 200,
                "base_condition": "lpc_contralign_ath02",
            },
            {
                "mode": "baseline",
                "budget": 300,
                "base_condition": "lpv_probcons_fpct05",
            },
        ],
        "fixed_parameters": {
            "beam": 100,
            "seed": 42,
            "omp_num_threads": 1,
            "refinement": 0,
            "viennarna": "2.7.2",
        },
        "quality_priority": ["mcc", "sci", "sps"],
        "descriptive_quality_thresholds": {"loose": 0.005, "strict": 0.001},
        "comparison_contract": {
            "candidate_budgets_vs_baseline500": [100, 200, 300],
            "baseline_budgets_vs_baseline500": [100, 200, 300],
            "same_budget_speed_gate": False,
            "tight_alignment_followup_eligible": False,
            "historical_max_iter_results_are_supplemental_only": True,
        },
        "interpretation": {
            "dd_iterations": "total iterations summed over progressive merges; max-iter is per merge",
            "dd_gap": "final dd_summary gap = best_ub - LB",
            "dd_violations": "dd_summary.violated at the last internal iterate of each merge",
            "validity": "metrics_validation and quality_complete are reported independently of DD disagreements",
        },
    }


def make_config(
    source: dict[str, Any],
    source_path: Path,
    binary: Path,
    results_root: Path,
    scope: str,
    mode: str,
    budget: int,
    dataset_ids: list[str],
    selected_conditions: list[str] | None = None,
) -> dict[str, Any]:
    config = copy.deepcopy(source)
    selected_conditions = selected_conditions or list(BASE_CONDITIONS)
    config["binary"] = str(binary.resolve())
    output_dir = results_root / scope / mode / f"iter-{budget}"
    if scope == "full-followup":
        output_dir /= selected_conditions[0]
    config["output_dir"] = str(
        output_dir.resolve()
    )
    source_map = source_conditions(source)
    config["conditions"] = [
        condition_for(source_map[name], mode, budget) for name in selected_conditions
    ]
    config["dataset_filter"] = {"ids": list(dataset_ids)}
    config["repetitions"] = 1
    config["random_seed"] = 42
    config["environment"] = {**config.get("environment", {}), "OMP_NUM_THREADS": "1"}
    # The source config has old experiment metadata that is not part of this
    # contract. Do not let it be mistaken for the DD-budget experiment.
    config.pop("comparison", None)
    config.pop("luna_relaxed", None)
    config["dd_budget"] = {
        "schema_version": 1,
        "scope": scope,
        "mode": mode,
        "mode_flag": MODE_FLAGS[mode],
        "complexity_note": MODE_COMPLEXITY[mode],
        "budget": budget,
        "budget_is_per_progressive_merge": True,
        "source_config": str(source_path.resolve()),
        "source_config_sha256": sha256_file(source_path),
        "binary_sha256": sha256_file(binary),
        "dataset_count": len(dataset_ids),
        "base_conditions": list(selected_conditions),
        "condition_count": len(selected_conditions),
        "expected_runs": len(dataset_ids) * len(selected_conditions),
        "post_pilot_condition_specific_followup": scope == "full-followup",
        "fixed_parameters": {
            "beam": 100,
            "seed": 42,
            "omp_num_threads": 1,
            "refinement": 0,
            "viennarna": "2.7.2",
        },
        "legacy_metadata_removed": True,
    }
    return config


def generated_specs(
    args: argparse.Namespace,
) -> list[tuple[str, str, int, list[str], list[str] | None]]:
    source = normalized_source(SOURCE_CONFIG)
    datasets = all_datasets(source)
    all_ids = [str(item["id"]) for item in datasets]
    if args.scope == "pilot":
        return [
            ("pilot", mode, budget, list(PILOT_DATASET_IDS), None)
            for mode in MODES
            for budget in BUDGETS
        ]
    if args.scope == "full-baseline":
        return [
            ("full-baseline", "baseline", budget, all_ids, None)
            for budget in (100, 500)
        ]
    if args.scope == "full-followup":
        if args.mode not in ("baseline", "projected-norm"):
            raise ValueError("full-followup requires baseline or projected-norm")
        if args.budget not in (200, 300):
            raise ValueError("full-followup budget must be 200 or 300")
        if args.base_condition not in BASE_CONDITIONS:
            raise ValueError(f"unknown base condition: {args.base_condition}")
        return [(
            "full-followup",
            args.mode,
            args.budget,
            all_ids,
            [args.base_condition],
        )]
    if args.mode not in ("projected-norm", "column-bound"):
        raise ValueError("full-candidate requires projected-norm or column-bound")
    if args.budget not in (100, 200, 300):
        raise ValueError("full-candidate budget must be 100, 200, or 300")
    return [("full-candidate", args.mode, args.budget, all_ids, None)]


def generate(args: argparse.Namespace) -> int:
    source = normalized_source(SOURCE_CONFIG)
    binary = absolute(args.binary, ROOT)
    if not binary.is_file():
        raise FileNotFoundError(binary)
    binary_sha = sha256_file(binary)
    if binary_sha != EXPECTED_BINARY_SHA256:
        raise RuntimeError(f"binary SHA256 mismatch: {binary_sha} != {EXPECTED_BINARY_SHA256}")
    results_root = absolute(args.results_root, ROOT)
    config_dir = absolute(args.config_dir, ROOT)
    plan_value = plan(binary)
    write_json(config_dir / "dd-budget-plan.json", plan_value, args.force)
    write_json(results_root / "plan.json", plan_value, args.force)
    for scope, mode, budget, dataset_ids, selected_conditions in generated_specs(args):
        config = make_config(
            source,
            SOURCE_CONFIG,
            binary,
            results_root,
            scope,
            mode,
            budget,
            dataset_ids,
            selected_conditions,
        )
        if scope == "pilot":
            name = f"dd-budget-pilot-{mode}-iter{budget}.json"
        elif scope == "full-baseline":
            name = f"dd-budget-full-baseline-iter{budget}.json"
        elif scope == "full-followup":
            name = (
                f"dd-budget-full-followup-{mode}-iter{budget}-"
                f"{selected_conditions[0]}.json"
            )
        else:
            name = f"dd-budget-full-candidate-{mode}-iter{budget}.json"
        path = config_dir / name
        write_json(path, config, args.force)
        print(
            f"{path} runs={len(dataset_ids) * len(config['conditions'])}"
        )
    print(f"plan={results_root / 'plan.json'}")
    return 0


def numeric(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def load_records(result_root: Path) -> list[dict[str, Any]]:
    """Read generated result trees using the existing trace_metrics helper."""
    if str(BENCHMARKS) not in sys.path:
        sys.path.insert(0, str(BENCHMARKS))
    from benchmark_luna_optimization import trace_metrics

    records: list[dict[str, Any]] = []
    for manifest_path in sorted(result_root.rglob("manifest.json")):
        manifest = read_json(manifest_path)
        config = manifest.get("config", {})
        metadata = config.get("dd_budget")
        if not isinstance(metadata, dict):
            continue
        condition_map = {
            str(condition["id"]): str(condition.get("base_condition", condition["id"]))
            for condition in config.get("conditions", [])
        }
        output_dir = manifest_path.parent
        for result_path in sorted(output_dir.glob("runs/*/*/rep-*/result.json")):
            result = read_json(result_path)
            metrics_path = result_path.parent / "metrics.jsonl"
            events = []
            if metrics_path.is_file():
                for line in metrics_path.read_text(encoding="utf-8").splitlines():
                    try:
                        events.append(json.loads(line))
                    except json.JSONDecodeError:
                        pass
            trace = trace_metrics(events, result)
            quality = result.get("quality") or {}
            condition = str(result.get("condition", ""))
            row: dict[str, Any] = {
                "scope": metadata.get("scope"),
                "mode": metadata.get("mode"),
                "budget": int(metadata["budget"]),
                "condition": condition_map.get(condition, condition),
                "dataset": str(result.get("dataset", "")),
                "repetition": int(result.get("repetition", 1)),
                "status": result.get("status"),
                "metrics_valid": bool(result.get("metrics_validation", {}).get("valid")),
                "quality_complete": bool(result.get("quality_complete")),
                "score_error": result.get("score_error"),
            }
            for field in QUALITY_FIELDS:
                row[field] = numeric(quality.get(field))
            for field in SUMMARY_FIELDS:
                if field in QUALITY_FIELDS:
                    continue
                row[field] = numeric(trace.get(field))
            row["dd_stop_reasons"] = json.dumps(trace.get("dd_stop_reasons", {}), sort_keys=True)
            records.append(row)
    return records


def mean_field(rows: list[dict[str, Any]], field: str) -> float | None:
    values = [numeric(row.get(field)) for row in rows]
    values = [value for value in values if value is not None]
    return fmean(values) if values else None


def median_field(rows: list[dict[str, Any]], field: str) -> float | None:
    values = [numeric(row.get(field)) for row in rows]
    values = [value for value in values if value is not None]
    return median(values) if values else None


def side_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "n": len(rows),
        "status_counts": dict(Counter(str(row.get("status")) for row in rows)),
        "metrics_valid": sum(bool(row.get("metrics_valid")) for row in rows),
        "quality_complete": sum(bool(row.get("quality_complete")) for row in rows),
        "means": {field: mean_field(rows, field) for field in SUMMARY_FIELDS},
    }


def quantile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    values = sorted(values)
    position = (len(values) - 1) * q
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    fraction = position - lower
    return values[lower] + (values[upper] - values[lower]) * fraction


def quality_thresholds(mean_deltas: dict[str, float | None]) -> dict[str, Any]:
    thresholds = {}
    for name, limit in (("loose", 0.005), ("strict", 0.001)):
        thresholds[name] = {
            "limit": limit,
            "absolute_mean_deltas": {
                field: None if mean_deltas[field] is None else abs(mean_deltas[field])
                for field in QUALITY_FIELDS
            },
            "all_within": all(
                mean_deltas[field] is not None
                and abs(mean_deltas[field]) <= limit
                for field in QUALITY_FIELDS
            ),
            "near_or_better": all(
                mean_deltas[field] is not None and mean_deltas[field] >= -limit
                for field in QUALITY_FIELDS
            ),
        }
    return thresholds


def paired_summary(
    candidate: list[dict[str, Any]], reference: list[dict[str, Any]], label: str
) -> dict[str, Any]:
    key = lambda row: (row["condition"], row["dataset"], row["repetition"])
    left = {key(row): row for row in candidate}
    candidate_conditions = {str(row["condition"]) for row in candidate}
    # A condition-specific follow-up may be paired with the all-four-condition
    # baseline500 tree. Filter only by condition, not by dataset: missing
    # datasets must remain visible as missing pairs.
    matched_reference = [
        row for row in reference if str(row["condition"]) in candidate_conditions
    ]
    right = {key(row): row for row in matched_reference}
    common = sorted(set(left) & set(right))
    pairs = []
    for item in common:
        a, b = left[item], right[item]
        row = {"condition": item[0], "dataset": item[1], "repetition": item[2]}
        for field in SUMMARY_FIELDS:
            av, bv = numeric(a.get(field)), numeric(b.get(field))
            row[f"candidate_{field}"] = av
            row[f"reference_{field}"] = bv
            row[f"delta_{field}"] = None if av is None or bv is None else av - bv
        pairs.append(row)

    distributions: dict[str, Any] = {}
    worst_cases: dict[str, Any] = {}
    for field in SUMMARY_FIELDS:
        delta_rows = [
            (numeric(row[f"delta_{field}"]), row)
            for row in pairs
            if numeric(row[f"delta_{field}"]) is not None
        ]
        deltas = [value for value, _row in delta_rows]
        distributions[field] = {
            "n": len(deltas),
            "mean_delta": fmean(deltas) if deltas else None,
            "median_delta": median(deltas) if deltas else None,
            "min_delta": min(deltas) if deltas else None,
            "p10_delta": quantile(deltas, 0.10),
            "p90_delta": quantile(deltas, 0.90),
            "max_delta": max(deltas) if deltas else None,
            "negative_count": sum(value < 0 for value in deltas),
            "positive_count": sum(value > 0 for value in deltas),
            "zero_count": sum(value == 0 for value in deltas),
        }
        if deltas:
            worst = min(delta_rows, key=lambda item: item[0])[1]
            best = max(delta_rows, key=lambda item: item[0])[1]
            worst_cases[field] = {
                "min_delta": {
                    "condition": worst["condition"],
                    "dataset": worst["dataset"],
                    "repetition": worst["repetition"],
                    "delta": worst[f"delta_{field}"],
                },
                "max_delta": {
                    "condition": best["condition"],
                    "dataset": best["dataset"],
                    "repetition": best["repetition"],
                    "delta": best[f"delta_{field}"],
                },
            }

    quality_deltas = {
        field: distributions[field]["mean_delta"] for field in QUALITY_FIELDS
    }
    thresholds = quality_thresholds(quality_deltas)
    per_condition = []
    for condition in BASE_CONDITIONS:
        condition_pairs = [row for row in pairs if row["condition"] == condition]
        condition_deltas = {
            field: (
                fmean(
                    numeric(row[f"delta_{field}"])
                    for row in condition_pairs
                    if numeric(row[f"delta_{field}"]) is not None
                )
                if any(numeric(row[f"delta_{field}"]) is not None for row in condition_pairs)
                else None
            )
            for field in QUALITY_FIELDS
        }
        per_condition.append(
            {
                "condition": condition,
                "candidate_means": {
                    field: (
                        fmean(
                            numeric(row[f"candidate_{field}"])
                            for row in condition_pairs
                            if numeric(row[f"candidate_{field}"]) is not None
                        )
                        if any(
                            numeric(row[f"candidate_{field}"]) is not None
                            for row in condition_pairs
                        )
                        else None
                    )
                    for field in SUMMARY_FIELDS
                },
                "reference_means": {
                    field: (
                        fmean(
                            numeric(row[f"reference_{field}"])
                            for row in condition_pairs
                            if numeric(row[f"reference_{field}"]) is not None
                        )
                        if any(
                            numeric(row[f"reference_{field}"]) is not None
                            for row in condition_pairs
                        )
                        else None
                    )
                    for field in SUMMARY_FIELDS
                },
                "mean_deltas": condition_deltas,
                "thresholds": quality_thresholds(condition_deltas),
                "paired_count": len(condition_pairs),
            }
        )
    return {
        "label": label,
        "candidate_count": len(candidate),
        "reference_input_count": len(reference),
        "reference_count": len(matched_reference),
        "paired_count": len(pairs),
        "missing_candidate": len(set(right) - set(left)),
        "missing_reference": len(set(left) - set(right)),
        "candidate": side_summary(candidate),
        # Report the condition-filtered reference side that was eligible for
        # pairing.  Keep reference_input_count above so an unexpected extra
        # condition remains visible without contaminating means/counts.
        "reference": side_summary(matched_reference),
        "mean_deltas": quality_deltas,
        "thresholds": thresholds,
        "per_condition": per_condition,
        "distributions": distributions,
        "worst_cases": worst_cases,
        "pairs": pairs,
    }


def aggregate_by_condition(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, int, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[
            (
                str(row["scope"]),
                int(row["budget"]),
                str(row["mode"]),
                str(row["condition"]),
            )
        ].append(row)
    output = []
    for (scope, budget, mode, condition), items in sorted(groups.items()):
        output.append(
            {
                "scope": scope,
                "mode": mode,
                "budget": budget,
                "condition": condition,
                **side_summary(items),
            }
        )
    return output


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize(args: argparse.Namespace) -> int:
    result_root = absolute(args.results_root, ROOT)
    rows = load_records(result_root)
    if not rows:
        raise RuntimeError(f"no dd-budget result records found under {result_root}")
    by_side: dict[tuple[str, int, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_side[(str(row["scope"]), int(row["budget"]), str(row["mode"]))].append(row)
    comparisons: list[dict[str, Any]] = []

    def add_comparison(
        candidate_key: tuple[str, int, str],
        reference_key: tuple[str, int, str],
        label: str,
    ) -> None:
        candidate = by_side.get(candidate_key, [])
        reference = by_side.get(reference_key, [])
        if candidate and reference:
            comparisons.append(paired_summary(candidate, reference, label))

    for mode in MODES:
        for budget in (100, 200, 300):
            add_comparison(
                ("pilot", budget, mode),
                ("pilot", 500, "baseline"),
                f"pilot:{mode}:iter{budget}:vs-baseline500",
            )
    add_comparison(
        ("full-baseline", 100, "baseline"),
        ("full-baseline", 500, "baseline"),
        "full:baseline:iter100:vs-baseline500",
    )
    for (scope, budget, mode), _items in sorted(by_side.items()):
        if scope == "full-candidate":
            add_comparison(
                (scope, budget, mode),
                ("full-baseline", 500, "baseline"),
                f"full:{mode}:iter{budget}:vs-baseline500",
            )
        if scope == "full-followup":
            configured_condition = sorted(
                {
                    str(row["condition"])
                    for row in by_side[(scope, budget, mode)]
                }
            )
            condition_label = configured_condition[0] if configured_condition else "unknown"
            add_comparison(
                (scope, budget, mode),
                ("full-baseline", 500, "baseline"),
                f"full-followup:{mode}:iter{budget}:{condition_label}:vs-baseline500",
            )
            if mode != "baseline":
                add_comparison(
                    (scope, budget, mode),
                    ("full-followup", budget, "baseline"),
                    f"full-followup:{mode}:iter{budget}:{condition_label}:vs-same-budget-baseline",
                )

    all_ids = {
        str(item["id"])
        for item in all_datasets(normalized_source(SOURCE_CONFIG))
    }
    strata = {
        "pilot12": set(PILOT_DATASET_IDS),
        "nonpilot73": all_ids - set(PILOT_DATASET_IDS),
    }
    stratified_comparisons = []

    def add_stratified(
        candidate_key: tuple[str, int, str],
        reference_key: tuple[str, int, str],
        label: str,
    ) -> None:
        for stratum, ids in strata.items():
            candidate = [
                row for row in by_side.get(candidate_key, []) if row["dataset"] in ids
            ]
            reference = [
                row for row in by_side.get(reference_key, []) if row["dataset"] in ids
            ]
            if candidate and reference:
                result = paired_summary(candidate, reference, f"{label}:{stratum}")
                result["stratum"] = stratum
                stratified_comparisons.append(result)

    add_stratified(
        ("full-baseline", 100, "baseline"),
        ("full-baseline", 500, "baseline"),
        "full:baseline:iter100:vs-baseline500",
    )
    for (scope, budget, mode), _items in sorted(by_side.items()):
        if scope == "full-followup":
            add_stratified(
                (scope, budget, mode),
                ("full-baseline", 500, "baseline"),
                f"full-followup:{mode}:iter{budget}:vs-baseline500",
            )
            if mode != "baseline":
                add_stratified(
                    (scope, budget, mode),
                    ("full-followup", budget, "baseline"),
                    f"full-followup:{mode}:iter{budget}:vs-same-budget-baseline",
                )

    plan_path = result_root / "plan.json"
    report = {
        "schema_version": 1,
        "experiment": "murlet-dataset1-dd-budget",
        "result_root": str(result_root),
        "plan": read_json(plan_path) if plan_path.is_file() else None,
        "record_count": len(rows),
        "side_summaries": [
            {"scope": key[0], "budget": key[1], "mode": key[2], **side_summary(items)}
            for key, items in sorted(by_side.items())
        ],
        "per_condition_means": aggregate_by_condition(rows),
        "stratified_comparisons": stratified_comparisons,
        "comparison_contract": {
            "candidate_100_200_300_are_compared_to_baseline500": True,
            "baseline_100_200_300_are_compared_to_baseline500": True,
            "same_budget_speed_gate": False,
            "quality_thresholds_are_descriptive_not_equivalence_claims": True,
            "tight_alignment_is_reference_only": True,
            "full_followup_is_post_pilot_condition_specific": True,
            "pilot_cases_are_not_held_out_from_followup": True,
        },
        "comparisons": comparisons,
    }
    write_json(result_root / "dd_budget_report.json", report)
    write_csv(result_root / "dd_budget_per_condition.csv", report["per_condition_means"])
    pair_rows = []
    for comparison in comparisons:
        for pair in comparison["pairs"]:
            pair_rows.append({"comparison": comparison["label"], **pair})
    write_csv(result_root / "dd_budget_pairs.csv", pair_rows)
    print(f"records={len(rows)} comparisons={len(comparisons)}")
    print(f"report={result_root / 'dd_budget_report.json'}")
    return 0


def selection_from_report(report: dict[str, Any]) -> dict[str, Any]:
    expected = len(PILOT_DATASET_IDS) * len(BASE_CONDITIONS)
    candidates = []
    for comparison in report.get("comparisons", []):
        label = comparison.get("label", "")
        if not label.startswith("pilot:") or ":vs-baseline500" not in label:
            continue
        parts = label.split(":")
        if (
            len(parts) != 4
            or parts[0] != "pilot"
            or parts[1] not in ("projected-norm", "column-bound")
            or parts[3] != "vs-baseline500"
        ):
            continue
        budget = int(parts[2].removeprefix("iter"))
        if budget not in (100, 200, 300):
            continue
        loose = comparison.get("thresholds", {}).get("loose", {})
        candidate_side = comparison.get("candidate", {})
        reference_side = comparison.get("reference", {})
        complete_side = lambda side: (
            side.get("n") == expected
            and side.get("metrics_valid") == expected
            and side.get("quality_complete") == expected
            and side.get("status_counts") == {"success": expected}
        )
        per_condition = comparison.get("per_condition", [])
        per_condition_ok = (
            len(per_condition) == len(BASE_CONDITIONS)
            and all(
                item.get("paired_count") == len(PILOT_DATASET_IDS)
                and item.get("thresholds", {})
                .get("loose", {})
                .get("near_or_better")
                for item in per_condition
            )
        )
        eligible = (
            bool(loose.get("near_or_better"))
            and comparison.get("paired_count") == expected
            and comparison.get("missing_candidate") == 0
            and comparison.get("missing_reference") == 0
            and complete_side(candidate_side)
            and complete_side(reference_side)
            and per_condition_ok
            and all(
                comparison.get("distributions", {})
                .get(field, {})
                .get("n") == expected
                for field in QUALITY_FIELDS
            )
        )
        if eligible:
            means = comparison.get("mean_deltas", {})
            candidates.append(
                {
                    "mode": parts[1],
                    "budget": budget,
                    "sci_delta": means.get("sci"),
                    "mcc_delta": means.get("mcc"),
                    "sps_delta": means.get("sps"),
                    "comparison": label,
                }
            )
    candidates.sort(
        key=lambda item: (
            item["budget"],
            -(item["sci_delta"] if item["sci_delta"] is not None else -float("inf")),
            -(item["mcc_delta"] if item["mcc_delta"] is not None else -float("inf")),
            -(item["sps_delta"] if item["sps_delta"] is not None else -float("inf")),
        )
    )
    selection = {
        "schema_version": 1,
        "rule": "eligible strictly-linear candidate has complete 48-pair success/valid pilot data, every combo is no worse than -0.005 mean MCC/SCI/SPS versus baseline500, and all overall means meet the same one-sided bound; choose lowest budget, then SCI, MCC, SPS",
        "eligible": candidates,
        "selected": candidates[0] if candidates else None,
        "full_baseline_still_required": True,
        "tight_alignment_excluded": True,
    }
    return selection


def select_candidate(args: argparse.Namespace) -> int:
    result_root = absolute(args.results_root, ROOT)
    report = read_json(result_root / "dd_budget_report.json")
    selection = selection_from_report(report)
    write_json(result_root / "dd_budget_selection.json", selection)
    print(json.dumps(selection, indent=2, ensure_ascii=False))
    return 0


SERIAL_DATASET_IDS = (
    "murlet-dataset1-RF00169-2",
    "murlet-dataset1-RF00017-1",
)


def bulk_prediction_path(
    bulk_root: Path, base_condition: str, variant: str, dataset_id: str
) -> Path:
    """Return the matching rep-01 prediction from the completed bulk arm."""
    specs = {
        "baseline500": ("full-baseline", "baseline", 500),
        "baseline100": ("full-baseline", "baseline", 100),
        "baseline200": ("full-followup", "baseline", 200),
        "baseline300": ("full-followup", "baseline", 300),
        "projected-norm200": ("full-followup", "projected-norm", 200),
        "projected-norm300": ("full-followup", "projected-norm", 300),
    }
    try:
        scope, mode, budget = specs[variant]
    except KeyError as error:
        raise ValueError(f"no bulk prediction mapping for serial variant {variant}") from error
    path = bulk_root / scope / mode / f"iter-{budget}"
    if scope == "full-followup":
        path /= base_condition
    return (
        path
        / "runs"
        / f"{base_condition}_iter{budget}"
        / dataset_id
        / "rep-01"
        / "prediction.fa"
    )


def serial(args: argparse.Namespace) -> int:
    """Run a small fixed-CPU, interleaved timing arm with differing budgets."""
    root = absolute(args.results_root, ROOT)
    if root.exists():
        raise RuntimeError(f"refusing to reuse serial result root: {root}")
    binary = absolute(args.binary, ROOT)
    if sha256_file(binary) != EXPECTED_BINARY_SHA256:
        raise RuntimeError("serial binary SHA256 does not match the experiment binary")
    affinity = os.sched_getaffinity(0)
    if args.cpu not in affinity:
        raise RuntimeError(f"requested serial CPU {args.cpu} is not in affinity {sorted(affinity)}")
    taskset = shutil.which("taskset")
    time_binary = shutil.which("time")
    if taskset is None or time_binary is None:
        raise RuntimeError("serial timing requires taskset and GNU time")

    source = normalized_source(SOURCE_CONFIG)
    source_map = source_conditions(source)
    manifest_path = Path(source["dataset_manifests"][0])
    dataset_map = {
        str(item["id"]): item for item in all_datasets(source)
    }
    metrics_enabled = not args.metrics_off
    bulk_root = root.parent
    variants = {
        "lpc_contralign_ath02": [
            ("baseline500", "baseline", 500),
            ("projected-norm200", "projected-norm", 200),
            ("baseline200", "baseline", 200),
            ("baseline100", "baseline", 100),
        ],
        "lpv_probcons_fpct05": [
            ("baseline500", "baseline", 500),
            ("projected-norm300", "projected-norm", 300),
            ("baseline300", "baseline", 300),
            ("baseline100", "baseline", 100),
        ],
    }
    root.mkdir(parents=True)
    records: list[dict[str, Any]] = []
    for repetition in range(1, args.repetitions + 1):
        for case_index, dataset_id in enumerate(SERIAL_DATASET_IDS):
            dataset = dataset_map[dataset_id]
            input_path = (manifest_path.parent / dataset["input"]).resolve()
            for condition_index, (base_condition, arms) in enumerate(variants.items()):
                rotation = (repetition + case_index + condition_index) % len(arms)
                ordered_arms = arms[rotation:] + arms[:rotation]
                for order_index, (variant, mode, budget) in enumerate(ordered_arms, 1):
                    run_dir = root / "runs" / variant / base_condition / dataset_id / f"rep-{repetition:02d}"
                    run_dir.mkdir(parents=True, exist_ok=True)
                    metrics_path = run_dir / "metrics.jsonl"
                    prediction_path = run_dir / "prediction.fa"
                    resource_path = run_dir / "resource.json"
                    stderr_path = run_dir / "stderr.log"
                    condition = condition_for(source_map[base_condition], mode, budget)
                    command = [
                        str(binary),
                        *condition["args"],
                    ]
                    if metrics_enabled:
                        command.extend(("--metrics-jsonl", str(metrics_path)))
                    command.append(str(input_path))
                    time_format = (
                        '{"wall_seconds":%e,"user_seconds":%U,'
                        '"system_seconds":%S,"max_rss_kb":%M,"exit_status":%x}'
                    )
                    wrapped = [
                        taskset,
                        "--cpu-list",
                        str(args.cpu),
                        time_binary,
                        "-f",
                        time_format,
                        "-o",
                        str(resource_path),
                        "--",
                        *command,
                    ]
                    started = time.monotonic()
                    with prediction_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
                        completed = subprocess.run(
                            wrapped,
                            stdout=stdout,
                            stderr=stderr,
                            env={**os.environ, "OMP_NUM_THREADS": "1"},
                            timeout=args.timeout_seconds,
                            check=False,
                        )
                    elapsed = time.monotonic() - started
                    resource = {}
                    if resource_path.is_file():
                        try:
                            resource = json.loads(resource_path.read_text(encoding="utf-8"))
                        except json.JSONDecodeError:
                            resource = {}
                    prediction_sha = (
                        sha256_file(prediction_path) if prediction_path.is_file() else None
                    )
                    matching_bulk_path = bulk_prediction_path(
                        bulk_root, base_condition, variant, dataset_id
                    )
                    bulk_prediction_sha = (
                        sha256_file(matching_bulk_path)
                        if matching_bulk_path.is_file()
                        else None
                    )
                    records.append(
                        {
                            "repetition": repetition,
                            "condition": base_condition,
                            "dataset": dataset_id,
                            "variant": variant,
                            "mode": mode,
                            "budget": budget,
                            "order_index": order_index,
                            "cpu": args.cpu,
                            "return_code": completed.returncode,
                            "status": "success" if completed.returncode == 0 else "failed",
                            "wall_seconds": numeric(resource.get("wall_seconds")),
                            "user_seconds": numeric(resource.get("user_seconds")),
                            "system_seconds": numeric(resource.get("system_seconds")),
                            "cpu_seconds": (
                                (numeric(resource.get("user_seconds")) or 0.0)
                                + (numeric(resource.get("system_seconds")) or 0.0)
                            ),
                            "max_rss_kb": numeric(resource.get("max_rss_kb")),
                            "harness_elapsed_seconds": elapsed,
                            "metrics_enabled": metrics_enabled,
                            "prediction_sha256": prediction_sha,
                            "bulk_prediction_path": str(matching_bulk_path),
                            "bulk_prediction_sha256": bulk_prediction_sha,
                            "prediction_matches_bulk": (
                                prediction_sha is not None
                                and prediction_sha == bulk_prediction_sha
                            ),
                            "command": command,
                        }
                    )
    write_json(
        root / "serial_manifest.json",
        {
            "schema_version": 1,
            "binary": str(binary),
            "binary_sha256": sha256_file(binary),
            "cpu": args.cpu,
            "omp_num_threads": 1,
            "repetitions": args.repetitions,
            "cases": list(SERIAL_DATASET_IDS),
            "metrics_enabled": metrics_enabled,
            "mode_label": "metrics-on" if metrics_enabled else "metrics-off",
            "bulk_result_root": str(bulk_root),
            "case_selection": (
                "fixed short/long representatives selected by input alignment length; "
                "selected before serial timing, after pilot; no performance-based selection"
            ),
            "interleaving": "variant order rotates within each repetition/case/condition",
            "prediction_validation": (
                "every serial prediction SHA256 must match the corresponding bulk "
                "arm's rep-01 prediction for the same case"
            ),
            "records": len(records),
        },
    )
    write_json(root / "serial_runs.json", records)
    write_csv(root / "serial_runs.csv", records)

    by_key = defaultdict(dict)
    for row in records:
        by_key[(row["condition"], row["dataset"], row["repetition"])][row["variant"]] = row
    pair_rows = []
    for key, arms in sorted(by_key.items()):
        candidate_variant = (
            "projected-norm200"
            if key[0] == "lpc_contralign_ath02"
            else "projected-norm300"
        )
        for reference_variant in (
            "baseline500",
            "baseline200" if key[0] == "lpc_contralign_ath02" else "baseline300",
            "baseline100",
        ):
            candidate = arms.get(candidate_variant)
            reference = arms.get(reference_variant)
            if not candidate or not reference:
                continue
            pair_rows.append(
                {
                    "condition": key[0],
                    "dataset": key[1],
                    "repetition": key[2],
                    "candidate": candidate_variant,
                    "reference": reference_variant,
                    "candidate_wall_seconds": candidate["wall_seconds"],
                    "reference_wall_seconds": reference["wall_seconds"],
                    "wall_delta_seconds": (
                        candidate["wall_seconds"] - reference["wall_seconds"]
                        if candidate["wall_seconds"] is not None
                        and reference["wall_seconds"] is not None
                        else None
                    ),
                    "wall_speedup_reference_over_candidate": (
                        reference["wall_seconds"] / candidate["wall_seconds"]
                        if candidate["wall_seconds"]
                        and reference["wall_seconds"]
                        else None
                    ),
                    "candidate_cpu_seconds": candidate["cpu_seconds"],
                    "reference_cpu_seconds": reference["cpu_seconds"],
                    "cpu_delta_seconds": candidate["cpu_seconds"] - reference["cpu_seconds"],
                    "cpu_speedup_reference_over_candidate": (
                        reference["cpu_seconds"] / candidate["cpu_seconds"]
                        if candidate["cpu_seconds"] > 0
                        else None
                    ),
                    "candidate_status": candidate["status"],
                    "reference_status": reference["status"],
                    "candidate_prediction_matches_bulk": candidate[
                        "prediction_matches_bulk"
                    ],
                    "reference_prediction_matches_bulk": reference[
                        "prediction_matches_bulk"
                    ],
                    "prediction_hashes_match_bulk": (
                        candidate["prediction_matches_bulk"]
                        and reference["prediction_matches_bulk"]
                    ),
                }
            )
    write_csv(root / "serial_pairs.csv", pair_rows)

    case_grouped = defaultdict(list)
    for row in pair_rows:
        case_grouped[
            (row["condition"], row["dataset"], row["candidate"], row["reference"])
        ].append(row)
    case_summary = []
    for key, items in sorted(case_grouped.items()):
        case_summary.append(
            {
                "condition": key[0],
                "dataset": key[1],
                "candidate": key[2],
                "reference": key[3],
                "n": len(items),
                "expected_repetitions": args.repetitions,
                "median_wall_delta_seconds": median_field(items, "wall_delta_seconds"),
                "median_wall_speedup": median_field(
                    items, "wall_speedup_reference_over_candidate"
                ),
                "median_cpu_delta_seconds": median_field(items, "cpu_delta_seconds"),
                "median_cpu_speedup": median_field(
                    items, "cpu_speedup_reference_over_candidate"
                ),
                "all_success": all(
                    row["candidate_status"] == row["reference_status"] == "success"
                    for row in items
                ),
                "all_prediction_hashes_match_bulk": all(
                    row["prediction_hashes_match_bulk"] for row in items
                ),
            }
        )

    summary = []
    grouped = defaultdict(list)
    for row in pair_rows:
        grouped[(row["condition"], row["candidate"], row["reference"])].append(row)
    for key, items in sorted(grouped.items()):
        summary.append(
            {
                "condition": key[0],
                "candidate": key[1],
                "reference": key[2],
                "n": len(items),
                "case_count": len({row["condition"] + ":" + row["dataset"] for row in items}),
                "mean_wall_delta_seconds": mean_field(items, "wall_delta_seconds"),
                "median_wall_delta_seconds": median_field(items, "wall_delta_seconds"),
                "mean_wall_speedup": mean_field(
                    items, "wall_speedup_reference_over_candidate"
                ),
                "median_wall_speedup": median_field(
                    items, "wall_speedup_reference_over_candidate"
                ),
                "mean_cpu_delta_seconds": mean_field(items, "cpu_delta_seconds"),
                "median_cpu_delta_seconds": median_field(items, "cpu_delta_seconds"),
                "mean_cpu_speedup": mean_field(
                    items, "cpu_speedup_reference_over_candidate"
                ),
                "median_cpu_speedup": median_field(
                    items, "cpu_speedup_reference_over_candidate"
                ),
                "all_success": all(
                    row["candidate_status"] == row["reference_status"] == "success"
                    for row in items
                ),
                "all_prediction_hashes_match_bulk": all(
                    row["prediction_hashes_match_bulk"] for row in items
                ),
            }
        )
    summary_payload = {
        "metrics_enabled": metrics_enabled,
        "rows": summary,
        "aggregate_rows": summary,
        "case_rows": case_summary,
        "pair_count": len(pair_rows),
        "aggregate_definition": (
            "aggregate wall/cpu speedups are means and medians of paired case-by-"
            "repetition ratios; aggregate seconds pool runs and are not a single "
            "case runtime. Per-case rows report medians over the five repetitions."
        ),
    }
    write_json(root / "serial_summary.json", summary_payload)
    write_csv(root / "serial_summary.csv", summary)
    write_csv(root / "serial_case_summary.csv", case_summary)
    print(
        f"serial_records={len(records)} pairs={len(pair_rows)} "
        f"metrics_enabled={metrics_enabled} root={root}"
    )
    all_success = bool(summary) and all(row["all_success"] for row in summary)
    all_hashes_match = bool(records) and all(
        row["prediction_matches_bulk"] for row in records
    )
    return 0 if all_success and all_hashes_match else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    plan_parser = sub.add_parser("plan", help="print the predeclared matrix")
    plan_parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)

    def handle_plan(args: argparse.Namespace) -> int:
        print(json.dumps(plan(absolute(args.binary, ROOT)), indent=2, ensure_ascii=False))
        return 0

    plan_parser.set_defaults(handler=handle_plan)

    generate_parser = sub.add_parser("generate", help="write configs for existing run.py")
    generate_parser.add_argument(
        "--scope",
        choices=("pilot", "full-baseline", "full-followup", "full-candidate"),
        required=True,
    )
    generate_parser.add_argument("--mode", choices=MODES)
    generate_parser.add_argument("--budget", type=int)
    generate_parser.add_argument("--base-condition", choices=BASE_CONDITIONS)
    generate_parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    generate_parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    generate_parser.add_argument("--config-dir", type=Path, default=BENCHMARKS / "configs")
    generate_parser.add_argument("--force", action="store_true")
    generate_parser.set_defaults(handler=generate)

    summary_parser = sub.add_parser("summarize", help="summarize completed result trees")
    summary_parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    summary_parser.set_defaults(handler=summarize)

    select_parser = sub.add_parser("select", help="select a strictly-linear pilot candidate")
    select_parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    select_parser.set_defaults(handler=select_candidate)

    serial_parser = sub.add_parser(
        "serial", help="run the predeclared fixed-CPU interleaved timing arm"
    )
    serial_parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    serial_parser.add_argument(
        "--results-root",
        type=Path,
        default=DEFAULT_RESULTS_ROOT / "serial-cpu19",
    )
    serial_parser.add_argument("--repetitions", type=int, default=5)
    serial_parser.add_argument("--timeout-seconds", type=float, default=3600.0)
    serial_parser.add_argument("--cpu", type=int, default=19)
    serial_parser.add_argument(
        "--metrics-off",
        action="store_true",
        help="omit --metrics-jsonl from timed commands and record metrics_enabled=false",
    )
    serial_parser.set_defaults(handler=serial)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return int(args.handler(args))
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
