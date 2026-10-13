#!/usr/bin/env python3
"""Generate small DD acceleration evaluation configs.

Execution remains with the existing ``benchmarks/run.py`` and Luna comparison
helpers.  This file only selects the pre-registered cases/conditions and
materializes separate old/off/on configs; it never builds or edits production.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"
SOURCE_CONFIG = BENCHMARKS / "configs/luna-relaxed.json"
BASELINE_BINARY = Path("/tmp/dafs-dd-accel-BpeZy7/dafs-baseline")
DEFAULT_CANDIDATE_BINARY = Path("/tmp/dafs-dd-accel-BpeZy7/dafs-candidate")
BASELINE_SHA256 = "47ae4d266b60a00291b71aab57abdc2175957d84fa324fc1e2e902170107da1e"
DEFAULT_FLAG = "--dd-projected-norm"
DEFAULT_RESULTS = BENCHMARKS / "results-murlet-dataset1-dd-acceleration-arm64-v272"

BASE_CONDITIONS = (
    "cl_linear_ath02",
    "lpv_probcons_fpct05",
    "lpc_contralign_ath02",
    "lpv_linear_fpct05",
)
PILOT_CONDITIONS = tuple(f"{name}_iter1000" for name in BASE_CONDITIONS)
FULL_CONDITIONS = tuple(
    condition_id
    for name in BASE_CONDITIONS
    for condition_id in (f"{name}_iter100", f"{name}_iter1000")
)
# One deterministic case from each of twelve families, spanning short to long
# alignments; RF00169-2 is also the independent CPU19 serial-A/B anchor.
PILOT_DATASETS = (
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


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any, force: bool) -> None:
    if path.exists() and not force:
        raise RuntimeError(f"refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(value: str | Path, base: Path = ROOT) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def datasets(config: dict[str, Any], config_path: Path) -> list[dict[str, Any]]:
    result = [copy.deepcopy(item) for item in config.get("datasets", [])]
    for manifest_name in config.get("dataset_manifests", []):
        manifest_path = resolve(manifest_name, config_path.parent)
        result.extend(read_json(manifest_path)["datasets"])
    selected = config.get("dataset_filter") or {}
    if selected.get("collections"):
        collections = {str(item) for item in selected["collections"]}
        result = [item for item in result if str(item.get("collection")) in collections]
    if selected.get("ids"):
        ids = {str(item) for item in selected["ids"]}
        result = [item for item in result if str(item.get("id")) in ids]
    return result


def normalize_paths(config: dict[str, Any], source_path: Path) -> dict[str, Any]:
    config = copy.deepcopy(config)
    base = source_path.parent
    config["dataset_manifests"] = [
        str(resolve(item, base)) for item in config.get("dataset_manifests", [])
    ]
    if config.get("sci", {}).get("bin_dir"):
        config["sci"]["bin_dir"] = str(resolve(config["sci"]["bin_dir"], base))
    return config


def make_config(
    source: dict[str, Any],
    source_path: Path,
    *,
    scope: str,
    variant: str,
    binary: Path,
    output_dir: Path,
    dataset_ids: list[str],
    condition_ids: list[str],
    candidate_flag: str,
    flag_on: bool,
) -> dict[str, Any]:
    config = normalize_paths(source, source_path)
    source_conditions = {str(item["id"]): item for item in config["conditions"]}
    missing_conditions = [item for item in condition_ids if item not in source_conditions]
    if missing_conditions:
        raise ValueError(f"missing source conditions: {missing_conditions}")
    available_datasets = {str(item["id"]) for item in datasets(config, source_path)}
    missing_datasets = [item for item in dataset_ids if item not in available_datasets]
    if missing_datasets:
        raise ValueError(f"missing source datasets: {missing_datasets}")
    selected_conditions = []
    for condition_id in condition_ids:
        condition = copy.deepcopy(source_conditions[condition_id])
        arguments = [str(item) for item in condition.get("args", [])]
        arguments = [item for item in arguments if item != candidate_flag]
        if flag_on:
            arguments.append(candidate_flag)
        condition["args"] = arguments
        selected_conditions.append(condition)
    config["binary"] = str(binary.resolve())
    config["output_dir"] = str(output_dir.resolve())
    config["conditions"] = selected_conditions
    config["dataset_filter"] = {"ids": list(dataset_ids)}
    config["repetitions"] = 1
    config["random_seed"] = 42
    config.setdefault("environment", {})["OMP_NUM_THREADS"] = "1"
    config["dd_acceleration"] = {
        "schema_version": 1,
        "driver": str(Path(__file__).resolve()),
        "source_config": str(source_path.resolve()),
        "source_config_sha256": sha256_file(source_path),
        "scope": scope,
        "variant": variant,
        "candidate_flag": candidate_flag,
        "flag_on": flag_on,
        "candidate_output_change_is_allowed": flag_on,
        "candidate_only_difference": "condition args contain the flag only when flag_on=true",
        "beam": 100,
        "pilot_dataset_count": len(PILOT_DATASETS),
        "dataset_count": len(dataset_ids),
        "condition_count": len(condition_ids),
        "expected_runs": len(dataset_ids) * len(condition_ids),
        "baseline_binary_sha256": BASELINE_SHA256,
    }
    if binary.is_file():
        config["dd_acceleration"]["binary_sha256"] = sha256_file(binary)
    # Keep the source's old comparison metadata from being mistaken for this
    # experiment's contract.
    config["comparison"] = {
        "mode": candidate_flag.removeprefix("--"),
        "candidate_only_arg": candidate_flag,
        "flag_on": flag_on,
        "output_difference_is_nonfatal": flag_on,
        "accuracy_filters_and_stop_criteria_changed": False,
    }
    return config


def generate(args: argparse.Namespace) -> int:
    source_path = resolve(args.source_config)
    source = read_json(source_path)
    if not isinstance(source, dict):
        raise ValueError(f"source config is not an object: {source_path}")
    if not args.candidate_flag.startswith("--") or "=" in args.candidate_flag:
        raise ValueError("--candidate-flag must be a bare long option")
    baseline = resolve(args.baseline_binary)
    candidate = resolve(args.candidate_binary)
    if not baseline.is_file():
        raise FileNotFoundError(f"baseline binary not found: {baseline}")
    baseline_sha = sha256_file(baseline)
    if baseline_sha != BASELINE_SHA256:
        raise RuntimeError(f"baseline SHA256 mismatch: {baseline_sha} != {BASELINE_SHA256}")
    all_dataset_ids = [str(item["id"]) for item in datasets(source, source_path)]
    pilot_ids = list(args.pilot_dataset_ids or PILOT_DATASETS)
    results_root = resolve(args.results_root)
    config_dir = resolve(args.config_dir)
    specs = (
        ("pilot-baseline", "pilot", "baseline", baseline, False, pilot_ids, list(PILOT_CONDITIONS), "pilot/baseline"),
        ("pilot-off-control", "pilot", "off-control", candidate, False, pilot_ids, list(PILOT_CONDITIONS), "pilot/off-control"),
        ("pilot-candidate", "pilot", "candidate", candidate, True, pilot_ids, list(PILOT_CONDITIONS), "pilot/candidate"),
        ("full-baseline", "full", "baseline", baseline, False, all_dataset_ids, list(FULL_CONDITIONS), "full/baseline"),
        ("full-candidate", "full", "candidate", candidate, True, all_dataset_ids, list(FULL_CONDITIONS), "full/candidate"),
    )
    paths: dict[str, str] = {}
    for name, scope, variant, binary, flag_on, dataset_ids, condition_ids, result_suffix in specs:
        config = make_config(
            source,
            source_path,
            scope=scope,
            variant=variant,
            binary=binary,
            output_dir=results_root / result_suffix,
            dataset_ids=dataset_ids,
            condition_ids=condition_ids,
            candidate_flag=args.candidate_flag,
            flag_on=flag_on,
        )
        path = config_dir / f"{args.prefix}-{name}.json"
        write_json(path, config, args.force)
        paths[name] = str(path)
        print(path)
    plan = {
        "schema_version": 1,
        "source_config": str(source_path),
        "source_config_sha256": sha256_file(source_path),
        "baseline_binary": str(baseline),
        "baseline_binary_sha256": baseline_sha,
        "candidate_binary": str(candidate),
        "candidate_binary_sha256": sha256_file(candidate) if candidate.is_file() else None,
        "candidate_flag": args.candidate_flag,
        "pilot_dataset_ids": pilot_ids,
        "pilot_expected_runs_per_tree": len(pilot_ids) * len(PILOT_CONDITIONS),
        "full_expected_runs_per_tree": len(all_dataset_ids) * len(FULL_CONDITIONS),
        "configs": paths,
        "execution": "generation only; wait for frozen candidate binary before run.py",
    }
    write_json(config_dir / f"{args.prefix}-plan.json", plan, args.force)
    print(config_dir / f"{args.prefix}-plan.json")
    print(
        f"prepared pilot={plan['pilot_expected_runs_per_tree']} "
        f"full={plan['full_expected_runs_per_tree']} runs/tree; no build or run executed"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("generate", nargs="?", default="generate", help=argparse.SUPPRESS)
    parser.add_argument("--source-config", type=Path, default=SOURCE_CONFIG)
    parser.add_argument("--config-dir", type=Path, default=BENCHMARKS / "configs")
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--baseline-binary", type=Path, default=BASELINE_BINARY)
    parser.add_argument("--candidate-binary", type=Path, default=DEFAULT_CANDIDATE_BINARY)
    parser.add_argument("--candidate-flag", default=DEFAULT_FLAG)
    parser.add_argument("--pilot-dataset-ids", nargs="+")
    parser.add_argument("--prefix", default="luna-dd-acceleration")
    parser.add_argument("--force", action="store_true")
    return parser


if __name__ == "__main__":
    try:
        raise SystemExit(generate(build_parser().parse_args()))
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        raise SystemExit(f"ERROR: {error}")
