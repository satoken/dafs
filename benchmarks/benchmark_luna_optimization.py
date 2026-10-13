#!/usr/bin/env python3
"""Generate and compare the reproducible Luna optimization benchmark.

This driver never builds or executes DAFS.  ``generate`` writes two configs
for the existing ``benchmarks/run.py`` runner; ``compare`` is intentionally
kept in this file so a result tree can be checked without rerunning it.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import re
import shlex
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import fmean
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = Path(__file__).resolve().parent
SOURCE_CONFIG = BENCHMARKS / "configs/murlet-dataset1-buffer-reuse-arm64-v272.json"
DEFAULT_BASELINE_BINARY = Path("/tmp/dafs-before-luna-20260919")
DEFAULT_OPTIMIZED_BINARY = ROOT / "build-arm64/src/dafs"
DEFAULT_HISTORICAL_RESULTS = (
    BENCHMARKS / "results-murlet-dataset1-buffer-reuse-arm64-v272"
)
MODES = ("decoder-only", "tight-alignment")
SUMMARY_TIME_KEYS = {"seconds", "dd_seconds"}
TRACE_RUN_START_IGNORED_KEYS = {
    "input",
    "binary",
    "binary_path",
    "binary_sha256",
    "binary_version",
    "executable",
    "executable_path",
}
METRIC_MEAN_FIELDS = (
    "mcc",
    "sci",
    "sps",
    "dd_iterations",
    "dd_wall_seconds",
    "wall_seconds",
    "user_seconds",
    "system_seconds",
    "cpu_seconds",
    "max_rss_kb",
    "dd_merge_count",
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

try:
    # Reuse the runner's metrics validation when comparison is requested.
    if str(BENCHMARKS) not in sys.path:
        sys.path.insert(0, str(BENCHMARKS))
    from run import load_metrics as runner_load_metrics
except ImportError:  # pragma: no cover - generate does not need the import
    runner_load_metrics = None


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def absolute_path(value: str | Path, base: Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def cli_path(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
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


def normalize_config_paths(config: dict[str, Any], source_path: Path) -> dict[str, Any]:
    """Make source-relative paths independent of the generated config location."""
    normalized = copy.deepcopy(config)
    source_base = source_path.parent
    if normalized.get("binary"):
        normalized["binary"] = str(absolute_path(normalized["binary"], source_base))
    if normalized.get("output_dir"):
        normalized["output_dir"] = str(
            absolute_path(normalized["output_dir"], source_base)
        )
    normalized["dataset_manifests"] = [
        str(absolute_path(value, source_base))
        for value in normalized.get("dataset_manifests", [])
    ]
    inline_datasets = normalized.get("datasets")
    if inline_datasets:
        for dataset in inline_datasets:
            if dataset.get("input"):
                dataset["input"] = str(absolute_path(dataset["input"], source_base))
            if dataset.get("reference"):
                dataset["reference"] = str(
                    absolute_path(dataset["reference"], source_base)
                )
    sci = normalized.get("sci")
    if isinstance(sci, dict):
        if sci.get("bin_dir"):
            sci["bin_dir"] = str(absolute_path(sci["bin_dir"], source_base))
        for key in ("rnaalifold", "rnafold"):
            if sci.get(key):
                sci[key] = str(absolute_path(sci[key], source_base))
    return normalized


def datasets_for_config(config: dict[str, Any], config_path: Path) -> list[dict[str, Any]]:
    """Load datasets using the same manifest/filter semantics as run.py."""
    datasets = [copy.deepcopy(item) for item in (config.get("datasets") or [])]
    config_base = config_path.parent
    for manifest_name in config.get("dataset_manifests", []) or []:
        manifest_path = absolute_path(manifest_name, config_base)
        manifest = read_json(manifest_path)
        datasets.extend(copy.deepcopy(manifest["datasets"]))
    dataset_filter = config.get("dataset_filter") or {}
    collections = dataset_filter.get("collections")
    if collections:
        allowed = {str(value) for value in collections}
        datasets = [
            item for item in datasets if str(item.get("collection")) in allowed
        ]
    ids = dataset_filter.get("ids")
    if ids:
        allowed = {str(value) for value in ids}
        datasets = [item for item in datasets if str(item.get("id")) in allowed]
    return datasets


def unique_dataset_ids(config: dict[str, Any], config_path: Path) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for dataset in datasets_for_config(config, config_path):
        dataset_id = str(dataset["id"])
        if dataset_id not in seen:
            seen.add(dataset_id)
            result.append(dataset_id)
    return result


def condition_max_iter(condition: dict[str, Any]) -> int | None:
    value = condition.get("max_iter")
    if value is not None:
        return int(value)
    arguments = list(map(str, condition.get("args", [])))
    for index, argument in enumerate(arguments):
        match = re.fullmatch(r"--max-iter(?:=(\d+))?", argument)
        if match:
            if match.group(1):
                return int(match.group(1))
            if index + 1 < len(arguments) and arguments[index + 1].isdigit():
                return int(arguments[index + 1])
    return None


def select_experiment(
    source: dict[str, Any],
    source_path: Path,
    max_iters: Iterable[int] | None = None,
    dataset_suffix: str | None = None,
    limit: int | None = None,
) -> dict[str, Any]:
    config = copy.deepcopy(source)
    if max_iters is not None:
        wanted = {int(value) for value in max_iters}
        config["conditions"] = [
            condition
            for condition in config.get("conditions", [])
            if condition_max_iter(condition) in wanted
        ]
        if not config["conditions"]:
            raise ValueError(f"no conditions selected for max-iter={sorted(wanted)}")

    if dataset_suffix is not None or limit is not None:
        dataset_ids = unique_dataset_ids(config, source_path)
        if dataset_suffix is not None:
            dataset_ids = [
                dataset_id for dataset_id in dataset_ids if dataset_id.endswith(dataset_suffix)
            ]
        if limit is not None:
            if limit <= 0:
                raise ValueError("--limit must be positive")
            dataset_ids = dataset_ids[:limit]
        if not dataset_ids:
            raise ValueError("dataset selection is empty")
        dataset_filter = copy.deepcopy(config.get("dataset_filter") or {})
        dataset_filter["ids"] = dataset_ids
        config["dataset_filter"] = dataset_filter
    return normalize_config_paths(config, source_path)


def add_tight_alignment_flag(config: dict[str, Any], mode: str) -> None:
    if mode != "tight-alignment":
        return
    for condition in config["conditions"]:
        args = [str(argument) for argument in condition.get("args", [])]
        if "--dd-tight-alignment" not in args:
            args.append("--dd-tight-alignment")
        condition["args"] = args


def set_sci_overrides(config: dict[str, Any], args: argparse.Namespace) -> None:
    if args.sci_bin_dir is None and args.sci_version is None:
        return
    sci = copy.deepcopy(config.get("sci") or {})
    if args.sci_bin_dir is not None:
        sci.pop("rnaalifold", None)
        sci.pop("rnafold", None)
        sci["bin_dir"] = str(cli_path(args.sci_bin_dir))
    if args.sci_version is not None:
        sci["version"] = args.sci_version
    config["sci"] = sci


def variant_config(
    base: dict[str, Any],
    source_path: Path,
    binary: Path,
    output_dir: Path,
    variant: str,
    mode: str,
    max_iters: list[int],
    dataset_suffix: str | None,
    limit: int | None,
) -> dict[str, Any]:
    config = normalize_config_paths(base, source_path)
    config["binary"] = str(binary)
    config["output_dir"] = str(output_dir)
    dataset_count = len(unique_dataset_ids(config, source_path))
    config["luna_optimization"] = {
        "driver": str(Path(__file__).resolve()),
        "source_config": str(source_path),
        "source_config_sha256": sha256_file(source_path),
        "variant": variant,
        "mode": mode,
        "dd_tight_alignment": mode == "tight-alignment",
        "max_iter_values": max_iters,
        "dataset_suffix": dataset_suffix,
        "dataset_limit": limit,
        "condition_count": len(config.get("conditions", [])),
        "dataset_count": dataset_count,
        "expected_runs": dataset_count * len(config.get("conditions", [])) * int(
            config.get("repetitions", 1)
        ),
    }
    if binary.is_file():
        config["luna_optimization"]["binary_sha256"] = sha256_file(binary)
    return config


def mode_options(parser: argparse.ArgumentParser) -> None:
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--mode", choices=MODES, default=None)
    group.add_argument(
        "--dd-tight-alignment",
        dest="mode",
        action="store_const",
        const="tight-alignment",
        help="append the experimental DD tight-alignment flag to every condition",
    )


def add_selection_options(
    parser: argparse.ArgumentParser,
    default_max_iters: list[int] | None,
) -> None:
    parser.add_argument(
        "--max-iter",
        dest="max_iters",
        type=int,
        nargs="+",
        choices=(100, 1000),
        default=default_max_iters,
        metavar="N",
        help="select max-iter conditions (default: 100 1000 for generate)",
    )
    parser.add_argument(
        "--dataset-suffix",
        help="keep dataset IDs ending with this suffix, e.g. --dataset-suffix=-0",
    )
    parser.add_argument("--limit", type=int, help="keep only the first selected datasets")


def generate(args: argparse.Namespace) -> int:
    source_path = cli_path(args.source_config)
    source = read_json(source_path)
    if not isinstance(source, dict):
        raise ValueError(f"source config is not an object: {source_path}")
    mode = args.mode or "decoder-only"
    max_iters = [int(value) for value in (args.max_iters or [100, 1000])]
    base = select_experiment(
        source,
        source_path,
        max_iters=max_iters,
        dataset_suffix=args.dataset_suffix,
        limit=args.limit,
    )
    add_tight_alignment_flag(base, mode)
    set_sci_overrides(base, args)

    config_dir = cli_path(args.config_dir)
    if args.results_root is None:
        result_name = (
            "results-murlet-dataset1-luna-decoder-arm64-v272"
            if mode == "decoder-only"
            else "results-murlet-dataset1-luna-tight-alignment-arm64-v272"
        )
        results_root = BENCHMARKS / result_name
    else:
        results_root = cli_path(args.results_root)
    baseline_results = cli_path(args.baseline_results) if args.baseline_results else results_root / "baseline"
    optimized_results = (
        cli_path(args.optimized_results) if args.optimized_results else results_root / "optimized"
    )
    baseline_binary = cli_path(args.baseline_binary)
    optimized_binary = cli_path(args.optimized_binary)
    prefix = args.prefix or ("luna-decoder" if mode == "decoder-only" else "luna-tight")

    baseline_config = variant_config(
        base,
        source_path,
        baseline_binary,
        baseline_results,
        "baseline",
        mode,
        max_iters,
        args.dataset_suffix,
        args.limit,
    )
    optimized_config = variant_config(
        base,
        source_path,
        optimized_binary,
        optimized_results,
        "optimized",
        mode,
        max_iters,
        args.dataset_suffix,
        args.limit,
    )
    baseline_path = config_dir / f"{prefix}-baseline.json"
    optimized_path = config_dir / f"{prefix}-optimized.json"
    write_json(baseline_path, baseline_config)
    write_json(optimized_path, optimized_config)

    dataset_count = int(optimized_config["luna_optimization"]["dataset_count"])
    condition_count = len(optimized_config["conditions"])
    expected_runs = int(optimized_config["luna_optimization"]["expected_runs"])
    print(
        f"generated mode={mode} conditions={condition_count} datasets={dataset_count} "
        f"runs_per_variant={expected_runs}"
    )
    print(f"baseline_config={baseline_path}")
    print(f"optimized_config={optimized_path}")
    print("run commands (not executed):")
    print(f"  python3 benchmarks/run.py --config {shlex.quote(str(baseline_path))}")
    print(f"  python3 benchmarks/run.py --config {shlex.quote(str(optimized_path))}")
    return 0


def key_text(key: tuple[str, str, int]) -> str:
    condition, dataset, repetition = key
    return f"{condition}/{dataset}/rep-{repetition:02d}"


def problem(side: str, key: tuple[str, str, int] | None, message: str) -> dict[str, Any]:
    return {
        "side": side,
        "key": key_text(key) if key is not None else None,
        "message": message,
    }


def discover_results(
    root: Path, side: str
) -> tuple[dict[tuple[str, str, int], dict[str, Any]], list[dict[str, Any]]]:
    found: dict[tuple[str, str, int], dict[str, Any]] = {}
    errors: list[dict[str, Any]] = []
    runs_root = root / "runs"
    if not runs_root.is_dir():
        errors.append(problem(side, None, f"missing result tree: {runs_root}"))
        return found, errors
    for path in sorted(runs_root.glob("*/*/rep-*/result.json")):
        try:
            relative = path.relative_to(runs_root)
            if len(relative.parts) != 4 or relative.parts[-1] != "result.json":
                raise ValueError(f"unexpected result path: {path}")
            repetition_match = re.fullmatch(r"rep-(\d+)", relative.parts[2])
            if repetition_match is None:
                raise ValueError(f"invalid repetition directory: {path}")
            key = (
                relative.parts[0],
                relative.parts[1],
                int(repetition_match.group(1)),
            )
            if key in found:
                errors.append(problem(side, key, f"duplicate result: {path}"))
                continue
            result = read_json(path)
            if not isinstance(result, dict):
                raise ValueError("result.json is not an object")
            found[key] = {"path": path, "run_dir": path.parent, "result": result}
        except (OSError, ValueError, json.JSONDecodeError) as error:
            errors.append(problem(side, None, f"{path}: {error}"))
    return found, errors


def event_stream(path: Path, errors: list[str]):
    """Yield parsed events without retaining a metrics file in memory."""
    if not path.is_file():
        errors.append(f"missing metrics stream: {path}")
        return
    try:
        with path.open(encoding="utf-8") as stream:
            for number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                try:
                    event = json.loads(line)
                    if not isinstance(event, dict):
                        raise ValueError("event is not an object")
                    yield event
                except (ValueError, json.JSONDecodeError) as error:
                    errors.append(f"{path}: line {number}: {error}")
    except OSError as error:
        errors.append(f"could not read metrics stream {path}: {error}")


def empty_trace_accumulator() -> dict[str, Any]:
    return {"dd_summaries": [], "last_iterations": {}, "run_summary": {}}


def consume_trace_event(accumulator: dict[str, Any], event: dict[str, Any]) -> None:
    event_name = event.get("event")
    if event_name == "dd_summary":
        accumulator["dd_summaries"].append(event)
    elif event_name == "dd_iteration" and event.get("merge_id") is not None:
        accumulator["last_iterations"][str(event["merge_id"])] = event
    elif event_name == "run_summary":
        accumulator["run_summary"] = event


def finish_trace_accumulator(
    accumulator: dict[str, Any], result: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    # Only the events needed by the requested aggregates survive; full event
    # equivalence is checked in the streaming loop below.
    compact_events = list(accumulator["dd_summaries"])
    compact_events.extend(
        event for _, event in sorted(accumulator["last_iterations"].items())
    )
    if accumulator["run_summary"]:
        compact_events.append(accumulator["run_summary"])
    summary = accumulator["run_summary"] or result.get("metrics_validation", {}).get(
        "run_summary", {}
    ) or {}
    return trace_metrics(compact_events, result), summary


def canonical_event(event: dict[str, Any]) -> dict[str, Any]:
    item = dict(event)
    item.pop("seconds", None)
    item.pop("dd_seconds", None)
    if item.get("event") == "run_start":
        for key in TRACE_RUN_START_IGNORED_KEYS:
            item.pop(key, None)
    return item


def scan_trace_pair(
    baseline_path: Path,
    optimized_path: Path,
    baseline_result: dict[str, Any],
    optimized_result: dict[str, Any],
    compare_traces: bool,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any], bool | None, dict[str, Any] | None, list[dict[str, Any]]]:
    baseline_errors: list[str] = []
    optimized_errors: list[str] = []
    baseline_stream = event_stream(baseline_path, baseline_errors)
    optimized_stream = event_stream(optimized_path, optimized_errors)
    baseline_accumulator = empty_trace_accumulator()
    optimized_accumulator = empty_trace_accumulator()
    sentinel = object()
    baseline_event = next(baseline_stream, sentinel)
    optimized_event = next(optimized_stream, sentinel)
    event_index = 0
    traces_equal = True if compare_traces else None
    first_trace_difference: dict[str, Any] | None = None
    while baseline_event is not sentinel or optimized_event is not sentinel:
        if baseline_event is not sentinel:
            consume_trace_event(baseline_accumulator, baseline_event)
        if optimized_event is not sentinel:
            consume_trace_event(optimized_accumulator, optimized_event)
        if compare_traces:
            if baseline_event is sentinel or optimized_event is sentinel:
                traces_equal = False
                if first_trace_difference is None:
                    first_trace_difference = {
                        "index": event_index,
                        "baseline_event": (
                            baseline_event.get("event")
                            if isinstance(baseline_event, dict)
                            else None
                        ),
                        "optimized_event": (
                            optimized_event.get("event")
                            if isinstance(optimized_event, dict)
                            else None
                        ),
                    }
            elif canonical_event(baseline_event) != canonical_event(optimized_event):
                traces_equal = False
                if first_trace_difference is None:
                    first_trace_difference = {
                        "index": event_index,
                        "baseline_event": baseline_event.get("event"),
                        "optimized_event": optimized_event.get("event"),
                    }
        event_index += 1
        baseline_event = next(baseline_stream, sentinel)
        optimized_event = next(optimized_stream, sentinel)

    errors: list[dict[str, Any]] = []
    errors.extend(problem("baseline", None, message) for message in baseline_errors)
    errors.extend(problem("optimized", None, message) for message in optimized_errors)
    baseline_metrics, baseline_summary = finish_trace_accumulator(
        baseline_accumulator, baseline_result
    )
    optimized_metrics, optimized_summary = finish_trace_accumulator(
        optimized_accumulator, optimized_result
    )
    if compare_traces and errors:
        traces_equal = None
    return (
        baseline_metrics,
        optimized_metrics,
        baseline_summary,
        optimized_summary,
        traces_equal,
        first_trace_difference,
        errors,
    )


def numeric(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def mean_or_none(values: Iterable[Any]) -> float | None:
    numbers = [value for item in values if (value := numeric(item)) is not None]
    return fmean(numbers) if numbers else None


def trace_metrics(events: list[dict[str, Any]], result: dict[str, Any]) -> dict[str, Any]:
    summaries = [event for event in events if event.get("event") == "dd_summary"]
    iterations_by_merge: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for event in events:
        if event.get("event") == "dd_iteration" and event.get("merge_id") is not None:
            iterations_by_merge[str(event["merge_id"])].append(event)
    last_iterations = [items[-1] for _, items in sorted(iterations_by_merge.items()) if items]
    run_summary = next(
        (event for event in reversed(events) if event.get("event") == "run_summary"),
        {},
    )
    if not run_summary:
        run_summary = result.get("metrics_validation", {}).get("run_summary", {}) or {}

    dd_iterations = numeric(run_summary.get("dd_iterations"))
    if dd_iterations is None:
        dd_iterations = sum(
            value
            for event in summaries
            if (value := numeric(event.get("iterations"))) is not None
        )
        if not dd_iterations and events:
            dd_iterations = float(sum(len(items) for items in iterations_by_merge.values()))

    violations = [numeric(event.get("violated")) for event in summaries]
    violations = [value for value in violations if value is not None]
    gaps = [numeric(event.get("gap")) for event in summaries]
    gaps = [value for value in gaps if value is not None]
    certified_violations = [
        numeric(event.get("certified_track_violated")) for event in last_iterations
    ]
    certified_violations = [value for value in certified_violations if value is not None]
    certified_gaps = [numeric(event.get("gap")) for event in last_iterations]
    certified_gaps = [value for value in certified_gaps if value is not None]
    stop_reasons = Counter(str(event.get("stop_reason", "<missing>")) for event in summaries)

    resource = result.get("resource", {}) or {}
    user_seconds = numeric(resource.get("user_seconds"))
    system_seconds = numeric(resource.get("system_seconds"))
    cpu_seconds = None
    if user_seconds is not None or system_seconds is not None:
        cpu_seconds = (user_seconds or 0.0) + (system_seconds or 0.0)
    quality = result.get("quality") or {}
    return {
        "mcc": numeric(quality.get("mcc")),
        "sci": numeric(quality.get("sci")),
        "sps": numeric(quality.get("sps")),
        "dd_iterations": dd_iterations,
        # DAFS emits this from steady_clock; it is wall time, not CPU time.
        "dd_wall_seconds": numeric(run_summary.get("dd_seconds")),
        "wall_seconds": numeric(resource.get("wall_seconds")),
        "user_seconds": user_seconds,
        "system_seconds": system_seconds,
        "cpu_seconds": cpu_seconds,
        "max_rss_kb": numeric(resource.get("max_rss_kb")),
        "dd_merge_count": len(summaries),
        "dd_violations_total": sum(violations) if violations else None,
        "dd_violations_mean_merge": mean_or_none(violations),
        "dd_violations_max_merge": max(violations) if violations else None,
        "dd_final_gap_mean": mean_or_none(gaps),
        "dd_final_gap_last": gaps[-1] if gaps else None,
        "dd_final_gap_max": max(gaps) if gaps else None,
        "certified_violations_total": (
            sum(certified_violations) if certified_violations else None
        ),
        "certified_violations_mean_merge": mean_or_none(certified_violations),
        "certified_final_gap_mean": mean_or_none(certified_gaps),
        "certified_final_gap_max": max(certified_gaps) if certified_gaps else None,
        "dd_stop_reasons": dict(sorted(stop_reasons.items())),
    }


def summary_without_time(summary: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in summary.items()
        if key not in SUMMARY_TIME_KEYS and not key.endswith("_seconds")
    }


def canonical_trace(events: list[dict[str, Any]]) -> list[dict[str, Any]]:
    canonical: list[dict[str, Any]] = []
    for event in events:
        item = dict(event)
        item.pop("seconds", None)
        item.pop("dd_seconds", None)
        if item.get("event") == "run_start":
            for key in TRACE_RUN_START_IGNORED_KEYS:
                item.pop(key, None)
        canonical.append(item)
    return canonical


def compare_pair(
    key: tuple[str, str, int],
    baseline: dict[str, Any],
    optimized: dict[str, Any],
    compare_traces: bool,
) -> dict[str, Any]:
    errors: list[dict[str, Any]] = []
    baseline_result = baseline["result"]
    optimized_result = optimized["result"]
    for side, info, result in (
        ("baseline", baseline, baseline_result),
        ("optimized", optimized, optimized_result),
    ):
        if result.get("status") != "success":
            errors.append(problem(side, key, f"status={result.get('status')}"))
        validation = result.get("metrics_validation")
        if not isinstance(validation, dict) or not validation.get("valid", False):
            errors.append(problem(side, key, "result metrics_validation.valid is false"))
        if result.get("score_error"):
            errors.append(problem(side, key, f"score_error={result['score_error']}"))

    (
        baseline_metrics,
        optimized_metrics,
        baseline_summary,
        optimized_summary,
        trace_equal,
        trace_difference,
        trace_errors,
    ) = scan_trace_pair(
        baseline["run_dir"] / "metrics.jsonl",
        optimized["run_dir"] / "metrics.jsonl",
        baseline_result,
        optimized_result,
        compare_traces,
    )
    errors.extend(
        {
            **error,
            "key": key_text(key),
        }
        for error in trace_errors
    )
    baseline_prediction = baseline["run_dir"] / "prediction.fa"
    optimized_prediction = optimized["run_dir"] / "prediction.fa"
    prediction_equal: bool | None = None
    if baseline_prediction.is_file() and optimized_prediction.is_file():
        prediction_equal = baseline_prediction.read_bytes() == optimized_prediction.read_bytes()
    else:
        if not baseline_prediction.is_file():
            errors.append(problem("baseline", key, f"missing prediction: {baseline_prediction}"))
        if not optimized_prediction.is_file():
            errors.append(problem("optimized", key, f"missing prediction: {optimized_prediction}"))

    summary_equal: bool | None = None
    if baseline_summary and optimized_summary:
        summary_equal = summary_without_time(baseline_summary) == summary_without_time(
            optimized_summary
        )
    else:
        errors.append(problem("comparison", key, "missing run_summary event"))

    return {
        "condition": key[0],
        "dataset": key[1],
        "repetition": key[2],
        "key": key_text(key),
        "prediction_bytes_equal": prediction_equal,
        "quality_dict_equal": baseline_result.get("quality") == optimized_result.get("quality"),
        "summary_non_time_equal": summary_equal,
        "trace_events_equal": trace_equal,
        "trace_first_difference": trace_difference,
        "baseline": baseline_metrics,
        "optimized": optimized_metrics,
        "errors": errors,
    }


def side_stats(records: list[dict[str, Any]], side: str) -> dict[str, Any]:
    metrics = [record[side] for record in records]
    means = {
        field: mean_or_none(metric.get(field) for metric in metrics)
        for field in METRIC_MEAN_FIELDS
    }
    stop_reasons: Counter[str] = Counter()
    for metric in metrics:
        stop_reasons.update(metric.get("dd_stop_reasons", {}))
    mode = None
    if stop_reasons:
        mode = max(stop_reasons.items(), key=lambda item: (item[1], item[0]))[0]
    return {
        "n": len(metrics),
        "means": means,
        "dd_stop_reasons": dict(sorted(stop_reasons.items())),
        "dd_stop_reason_mode": mode,
    }


def speedup(old: Any, new: Any) -> float | None:
    old_number, new_number = numeric(old), numeric(new)
    if old_number is None or new_number is None or new_number == 0:
        return None
    return old_number / new_number


def reduction_percent(old: Any, new: Any) -> float | None:
    old_number, new_number = numeric(old), numeric(new)
    if old_number is None or new_number is None or old_number == 0:
        return None
    return 100.0 * (old_number - new_number) / old_number


def change_percent(old: Any, new: Any) -> float | None:
    old_number, new_number = numeric(old), numeric(new)
    if old_number is None or new_number is None or old_number == 0:
        return None
    return 100.0 * (new_number - old_number) / old_number


def timing_comparison(baseline: dict[str, Any], optimized: dict[str, Any]) -> dict[str, Any]:
    old_means = baseline["means"]
    new_means = optimized["means"]
    result: dict[str, Any] = {}
    for field, label in (
        ("wall_seconds", "wall"),
        ("dd_wall_seconds", "dd_wall"),
        ("cpu_seconds", "cpu"),
    ):
        result[f"{label}_speedup"] = speedup(old_means.get(field), new_means.get(field))
        result[f"{label}_reduction_percent"] = reduction_percent(
            old_means.get(field), new_means.get(field)
        )
    result["rss_change_percent"] = change_percent(
        old_means.get("max_rss_kb"), new_means.get("max_rss_kb")
    )
    return result


def result_manifest_info(root: Path) -> dict[str, Any]:
    path = root / "manifest.json"
    if not path.is_file():
        return {"manifest_path": str(path), "present": False}
    try:
        manifest = read_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        return {"manifest_path": str(path), "present": False, "error": str(error)}
    return {
        "manifest_path": str(path),
        "present": True,
        "config_path": manifest.get("config_path"),
        "binary": manifest.get("binary"),
        "binary_sha256": manifest.get("binary_sha256"),
        "git_commit": manifest.get("git_commit"),
        "created_utc": manifest.get("created_utc"),
    }


def load_compare_config(
    explicit_path: Path | None, optimized_root: Path
) -> tuple[dict[str, Any], Path, list[str]]:
    warnings: list[str] = []
    manifest_path = optimized_root / "manifest.json"
    manifest: dict[str, Any] = {}
    if manifest_path.is_file():
        try:
            loaded = read_json(manifest_path)
            if isinstance(loaded, dict):
                manifest = loaded
        except (OSError, ValueError, json.JSONDecodeError) as error:
            warnings.append(f"could not read optimized manifest: {error}")

    if explicit_path is not None:
        config_path = cli_path(explicit_path)
        return read_json(config_path), config_path, warnings

    manifest_config_path = manifest.get("config_path")
    if manifest_config_path:
        config_path = cli_path(manifest_config_path)
        if config_path.is_file():
            return read_json(config_path), config_path, warnings
        warnings.append(f"manifest config_path does not exist: {config_path}")
    if isinstance(manifest.get("config"), dict):
        # The manifest writer preserves the original config path for resolving
        # any relative fields in this embedded copy.
        config_path = (
            cli_path(manifest_config_path)
            if manifest_config_path
            else SOURCE_CONFIG
        )
        return manifest["config"], config_path, warnings
    return read_json(SOURCE_CONFIG), SOURCE_CONFIG, warnings


def flatten_condition_report(report: dict[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "condition": report["condition"],
        "base_condition": report["base_condition"],
        "max_iter": report["max_iter"],
        "expected_runs": report["expected_runs"],
        "paired_runs": report["paired_runs"],
    }
    equality = report["equality"]
    row.update({f"{key}_equal_count": value for key, value in equality.items()})
    for side in ("baseline", "optimized"):
        stats = report[side]
        for field in METRIC_MEAN_FIELDS:
            row[f"{side}_{field}_mean"] = stats["means"].get(field)
        row[f"{side}_dd_stop_reasons"] = json.dumps(
            stats["dd_stop_reasons"], sort_keys=True, separators=(",", ":")
        )
        row[f"{side}_dd_stop_reason_mode"] = stats["dd_stop_reason_mode"]
    row.update(report["timing"])
    return row


def compare(args: argparse.Namespace) -> int:
    baseline_root = cli_path(args.baseline_results)
    optimized_root = cli_path(args.optimized_results)
    config, config_path, config_warnings = load_compare_config(args.config, optimized_root)
    selected = select_experiment(
        config,
        config_path,
        max_iters=args.max_iters,
        dataset_suffix=args.dataset_suffix,
        limit=args.limit,
    )
    mode = args.mode or selected.get("luna_optimization", {}).get("mode") or "decoder-only"
    if mode not in MODES:
        raise ValueError(f"unknown mode: {mode}")
    conditions = selected.get("conditions", [])
    dataset_ids = unique_dataset_ids(selected, config_path)
    repetitions = int(selected.get("repetitions", 1))
    expected_keys = {
        (str(condition["id"]), dataset_id, repetition)
        for dataset_id in dataset_ids
        for condition in conditions
        for repetition in range(1, repetitions + 1)
    }
    condition_metadata = {}
    for condition in conditions:
        condition_id = str(condition["id"])
        base_condition = condition.get("base_condition")
        if base_condition is None:
            base_condition = condition_id.rsplit("_iter", 1)[0]
        condition_metadata[condition_id] = {
            "base_condition": str(base_condition),
            "max_iter": condition_max_iter(condition),
        }

    baseline, baseline_errors = discover_results(baseline_root, "baseline")
    optimized, optimized_errors = discover_results(optimized_root, "optimized")
    errors = [*baseline_errors, *optimized_errors]
    missing_runs: list[dict[str, Any]] = []
    unexpected_runs: list[dict[str, Any]] = []
    for side, found in (("baseline", baseline), ("optimized", optimized)):
        for key in sorted(found.keys() - expected_keys):
            unexpected_runs.append({"side": side, "key": key_text(key)})
        for key in sorted(expected_keys - found.keys()):
            missing_runs.append({"side": side, "key": key_text(key)})

    pair_records: list[dict[str, Any]] = []
    mismatch_keys: dict[str, list[str]] = {
        "prediction_bytes": [],
        "quality_dict": [],
        "summary_non_time": [],
    }
    if args.compare_traces:
        mismatch_keys["trace_events"] = []
    for key in sorted(expected_keys & baseline.keys() & optimized.keys()):
        record = compare_pair(key, baseline[key], optimized[key], args.compare_traces)
        pair_records.append(record)
        errors.extend(record["errors"])
        if record["prediction_bytes_equal"] is False:
            mismatch_keys["prediction_bytes"].append(record["key"])
        if record["quality_dict_equal"] is False:
            mismatch_keys["quality_dict"].append(record["key"])
        if record["summary_non_time_equal"] is False:
            mismatch_keys["summary_non_time"].append(record["key"])
        if args.compare_traces and record["trace_events_equal"] is False:
            mismatch_keys["trace_events"].append(record["key"])

    by_condition: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in pair_records:
        by_condition[record["condition"]].append(record)
    condition_reports: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    for condition in conditions:
        condition_id = str(condition["id"])
        records = by_condition.get(condition_id, [])
        baseline_stats = side_stats(records, "baseline")
        optimized_stats = side_stats(records, "optimized")
        equality = {
            "prediction": sum(record["prediction_bytes_equal"] is True for record in records),
            "quality": sum(record["quality_dict_equal"] is True for record in records),
            "summary_non_time": sum(
                record["summary_non_time_equal"] is True for record in records
            ),
        }
        if args.compare_traces:
            equality["trace"] = sum(
                record["trace_events_equal"] is True for record in records
            )
        metadata = condition_metadata[condition_id]
        condition_report = {
            "condition": condition_id,
            "base_condition": metadata["base_condition"],
            "max_iter": metadata["max_iter"],
            "expected_runs": len(dataset_ids) * repetitions,
            "paired_runs": len(records),
            "equality": equality,
            "baseline": baseline_stats,
            "optimized": optimized_stats,
            "timing": timing_comparison(baseline_stats, optimized_stats),
        }
        condition_reports.append(condition_report)
        csv_rows.append(flatten_condition_report(condition_report))

    baseline_info = result_manifest_info(baseline_root)
    optimized_info = result_manifest_info(optimized_root)
    limitations = [
        "dd_wall_seconds is DAFS metrics run_summary.dd_seconds measured by steady_clock; it is wall time, not CPU time.",
        "Wall time, DD wall time, user+system CPU time, and RSS remain host measurements and do not isolate contention.",
    ]
    if args.baseline_kind == "historical":
        limitations.append(
            "The baseline result tree is historical rather than a fresh paired run; timing deltas are not a controlled same-window comparison."
        )
    else:
        limitations.append(
            "Fresh baseline/new pairing and any alternating schedule are external to this driver; the report matches keys but does not prove thermal or contention control."
        )
        limitations.append(
            "For a controlled timing subset, alternate matched baseline and optimized cases while using the same generated selection."
        )

    report_dir = (
        cli_path(args.report_dir)
        if args.report_dir is not None
        else optimized_root / "luna-optimization-report"
    )
    csv_path = report_dir / "luna_optimization_per_condition.csv"
    json_path = report_dir / "luna_optimization_report.json"
    report = {
        "schema_version": 1,
        "created_utc": now_utc(),
        "driver": str(Path(__file__).resolve()),
        "mode": mode,
        "compare_traces": bool(args.compare_traces),
        "baseline_kind": args.baseline_kind,
        "config": str(config_path),
        "experiment": {
            "condition_count": len(conditions),
            "dataset_count": len(dataset_ids),
            "repetitions": repetitions,
            "expected_runs": len(expected_keys),
            "conditions": [str(condition["id"]) for condition in conditions],
            "dataset_ids": dataset_ids,
        },
        "baseline": {"results": str(baseline_root), **baseline_info},
        "optimized": {"results": str(optimized_root), **optimized_info},
        "counts": {
            "expected_runs": len(expected_keys),
            "baseline_result_files": len(baseline),
            "optimized_result_files": len(optimized),
            "paired_runs": len(pair_records),
            "missing_runs": len(missing_runs),
            "unexpected_runs": len(unexpected_runs),
            "errors": len(errors),
            "prediction_bytes_equal": sum(
                record["prediction_bytes_equal"] is True for record in pair_records
            ),
            "quality_dict_equal": sum(
                record["quality_dict_equal"] is True for record in pair_records
            ),
            "summary_non_time_equal": sum(
                record["summary_non_time_equal"] is True for record in pair_records
            ),
            "trace_events_equal": (
                sum(record["trace_events_equal"] is True for record in pair_records)
                if args.compare_traces
                else None
            ),
        },
        "mismatches": mismatch_keys,
        "missing_runs": missing_runs,
        "unexpected_runs": unexpected_runs,
        "errors": errors,
        "warnings": config_warnings,
        "limitations": limitations,
        "per_condition": condition_reports,
        "per_run": pair_records,
    }
    report_dir.mkdir(parents=True, exist_ok=True)
    write_json(json_path, report)
    fieldnames = list(csv_rows[0]) if csv_rows else ["condition"]
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)

    decoder_mismatch = mode == "decoder-only" and any(mismatch_keys.values())
    # A larger result tree may intentionally contain runs outside a pilot's
    # selected subset.  Report those extras, but only missing/error evidence or
    # an actual decoder comparison mismatch makes the command fail.
    fatal = bool(errors or missing_runs or decoder_mismatch)
    print(
        f"comparison mode={mode} expected={len(expected_keys)} paired={len(pair_records)} "
        f"errors={len(errors)} missing={len(missing_runs)} unexpected={len(unexpected_runs)}"
    )
    print(f"csv_report={csv_path}")
    print(f"json_report={json_path}")
    if fatal:
        if decoder_mismatch:
            print("decoder mismatch: predictions/quality/summary (or requested traces) differ", file=sys.stderr)
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    generate_parser = subparsers.add_parser(
        "generate", help="write baseline and optimized configs; do not run them"
    )
    generate_parser.add_argument("--source-config", type=Path, default=SOURCE_CONFIG)
    generate_parser.add_argument(
        "--config-dir",
        "--out-dir",
        dest="config_dir",
        type=Path,
        default=BENCHMARKS / "configs",
        help="directory for generated JSON configs",
    )
    generate_parser.add_argument("--results-root", type=Path)
    generate_parser.add_argument("--baseline-results", type=Path)
    generate_parser.add_argument("--optimized-results", type=Path)
    generate_parser.add_argument(
        "--baseline-binary", type=Path, default=DEFAULT_BASELINE_BINARY
    )
    generate_parser.add_argument(
        "--optimized-binary", type=Path, default=DEFAULT_OPTIMIZED_BINARY
    )
    generate_parser.add_argument("--sci-bin-dir", type=Path)
    generate_parser.add_argument("--sci-version")
    generate_parser.add_argument("--prefix")
    mode_options(generate_parser)
    add_selection_options(generate_parser, [100, 1000])
    generate_parser.set_defaults(handler=generate)

    compare_parser = subparsers.add_parser(
        "compare", help="compare two completed result trees and write reports"
    )
    compare_parser.add_argument(
        "--baseline-results", type=Path, default=DEFAULT_HISTORICAL_RESULTS
    )
    compare_parser.add_argument(
        "--optimized-results",
        type=Path,
        default=BENCHMARKS / "results-murlet-dataset1-luna-decoder-arm64-v272/optimized",
    )
    compare_parser.add_argument("--config", type=Path)
    compare_parser.add_argument("--report-dir", type=Path)
    compare_parser.add_argument(
        "--baseline-kind",
        choices=("historical", "fresh"),
        default="historical",
        help="provenance note for timing interpretation (default: historical)",
    )
    compare_parser.add_argument(
        "--compare-traces",
        action="store_true",
        help="stream-compare all metrics events after documented normalization",
    )
    mode_options(compare_parser)
    add_selection_options(compare_parser, None)
    compare_parser.set_defaults(handler=compare)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return int(args.handler(args))
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
