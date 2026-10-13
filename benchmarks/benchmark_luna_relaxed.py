#!/usr/bin/env python3
"""Evaluation helpers for the approximate-output Luna relaxed experiment.

This driver owns only benchmark artifacts.  It never edits or rebuilds the
shared ``build-arm64`` tree.  Standalone decoder microbenchmarks are compiled
into an explicit result directory from either the frozen baseline source tree
or the candidate source tree.  Full DAFS runs are still performed by
``benchmarks/run.py`` using the checked-in JSON configs; ``compare`` treats
changed predictions and quality as measurements, not failures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"
SCALING_SOURCE = BENCHMARKS / "luna_round2_scaling.cpp"
DEFAULT_BASELINE_BINARY = Path("/tmp/dafs-luna-relaxed-X52WP9/dafs-baseline")
DEFAULT_BASELINE_SOURCE = Path("/tmp/dafs-luna-relaxed-X52WP9/src")
DEFAULT_CANDIDATE_SOURCE = ROOT / "src"
DEFAULT_BASELINE_RESULTS = (
    BENCHMARKS / "results-murlet-dataset1-luna-round2-arm64-v272"
)
DEFAULT_CANDIDATE_RESULTS = (
    BENCHMARKS / "results-murlet-dataset1-luna-relaxed-arm64-v272"
)
DEFAULT_CONFIG = BENCHMARKS / "configs/luna-relaxed.json"
DEFAULT_CPU = 19
DEFAULT_LENGTHS = (250, 500, 1000, 2000, 4000)
DEFAULT_VARIANTS = ("beam", "nuss-certified", "nuss-uncertified")
BASELINE_SHA256 = "75b6fa32805ec02cb3caf32186978cc6f6aed9bf0c51aec410a5347dda9a8531"
SERIAL_DATASET = "murlet-dataset1-RF00169-2"
SERIAL_CONDITION = "cl_linear_ath02_iter1000"

if str(BENCHMARKS) not in sys.path:
    sys.path.insert(0, str(BENCHMARKS))
from run import load_metrics as runner_load_metrics


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def resolve_path(value: str | Path, base: Path = ROOT) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def ensure_fresh(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise RuntimeError(f"refusing to overwrite non-empty path: {path}")
    path.mkdir(parents=True, exist_ok=True)


def source_manifest(
    source_root: Path, destination: Path, relative_paths: Iterable[str]
) -> str:
    """Hash only the decoder files used by the isolated microbenchmarks."""
    source_root = source_root.resolve()
    entries: list[str] = []
    for relative in sorted(set(relative_paths)):
        path = source_root / relative
        if not path.is_file():
            raise FileNotFoundError(f"microbenchmark source file not found: {path}")
        entries.append(f"{sha256_file(path)}  {relative}\n")
    destination.write_text("".join(entries), encoding="utf-8")
    return sha256_file(destination)


def host_context(cpu: int | None) -> dict[str, Any]:
    models: dict[str, str] = {}
    try:
        output = subprocess.check_output(
            ["lscpu", "-p=CPU,MODELNAME"], text=True, stderr=subprocess.DEVNULL
        )
        for line in output.splitlines():
            if not line or line.startswith("#") or "," not in line:
                continue
            number, model = line.split(",", 1)
            models[number] = model
    except (OSError, subprocess.CalledProcessError):
        pass
    return {
        "platform": platform.platform(),
        "uname": list(platform.uname()),
        "cpu_pin": cpu,
        "cpu_model_for_pin": models.get(str(cpu)) if cpu is not None else None,
        "cpu_models": models,
    }


def compiler_command(
    compiler: str, macro: str, source_root: Path, source: Path, output: Path
) -> list[str]:
    command = [compiler, "-O3", "-DNDEBUG", "-std=c++17", f"-D{macro}"]
    if macro == "LUNA_ROUND2_BEAM":
        command += ["-I", str(source_root / "linearalign"), "-I", str(source_root)]
    else:
        command += ["-I", str(source_root), "-include", str(source_root / "nussinov.h")]
    command += [str(SCALING_SOURCE), str(source), "-o", str(output)]
    return command


def compile_micro(args: argparse.Namespace) -> int:
    side = args.side
    source_root = resolve_path(
        args.source_root
        if args.source_root is not None
        else (DEFAULT_BASELINE_SOURCE if side == "baseline" else DEFAULT_CANDIDATE_SOURCE)
    )
    if not source_root.is_dir():
        raise FileNotFoundError(f"source root not found: {source_root}")
    root = resolve_path(args.root)
    ensure_fresh(root)
    executable_dir = root / "executables"
    executable_dir.mkdir()
    compiler = os.environ.get("CXX", "c++")
    specs = (
        ("beam", "LUNA_ROUND2_BEAM", source_root / "linearalign/BeamAlign.cpp"),
        ("nuss", "LUNA_ROUND2_NUSS", source_root / "nussinov.cpp"),
    )
    executables: dict[str, dict[str, Any]] = {}
    for name, macro, source in specs:
        if not source.is_file():
            raise FileNotFoundError(f"implementation source not found: {source}")
        output = executable_dir / f"{side}-{name}"
        command = compiler_command(compiler, macro, source_root, source, output)
        print("COMPILE", " ".join(command), flush=True)
        subprocess.run(command, cwd=ROOT, check=True)
        executables[name] = {
            "path": str(output),
            "sha256": sha256_file(output),
            "source": str(source),
            "source_sha256": sha256_file(source),
            "command": command,
        }
    manifest_path = root / "source-manifest.sha256"
    relevant_files = (
        "linearalign/BeamAlign.cpp",
        "linearalign/BeamAlign.h",
        "nussinov.cpp",
        "nussinov.h",
        "relaxed_bounds.h",
    )
    manifest_sha256 = source_manifest(source_root, manifest_path, relevant_files)
    write_json(
        root / "compile.json",
        {
            "schema_version": 1,
            "created_utc": utc_now(),
            "side": side,
            "compiler": compiler,
            "workspace": str(ROOT),
            "source_root": str(source_root),
            "source_manifest": str(manifest_path),
            "source_manifest_sha256": manifest_sha256,
            "source_manifest_files": list(relevant_files),
            "scaling_source": str(SCALING_SOURCE),
            "scaling_source_sha256": sha256_file(SCALING_SOURCE),
            "executables": executables,
            "host": host_context(None),
        },
    )
    print(root / "compile.json")
    return 0


def read_one_csv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise RuntimeError(f"expected one row in {path}, got {len(rows)}")
    return rows[0]


def run_command(command: list[str], cpu: int | None) -> list[str]:
    if cpu is None:
        return command
    if shutil.which("taskset") is None:
        raise RuntimeError("taskset is required for a pinned microbenchmark")
    return ["taskset", "-c", str(cpu), *command]


def read_micro_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def run_micro(args: argparse.Namespace) -> int:
    root = resolve_path(args.root)
    manifest = json.loads((root / "compile.json").read_text(encoding="utf-8"))
    side = str(manifest["side"])
    result_root = root / "results"
    ensure_fresh(result_root)
    lengths = tuple(int(value) for value in args.lengths)
    variants = tuple(args.variants)
    if not lengths or not variants:
        raise ValueError("microbenchmark selection is empty")
    if any(variant not in DEFAULT_VARIANTS for variant in variants):
        raise ValueError(f"unknown micro variant: {variants}")
    rows: list[dict[str, str]] = []
    invocation = 0
    for variant in variants:
        executable_name = "beam" if variant == "beam" else "nuss"
        executable = Path(manifest["executables"][executable_name]["path"])
        for length in lengths:
            invocation += 1
            row_path = result_root / f"{invocation:02d}-{side}-{variant}-{length}.csv"
            command = [
                str(executable),
                "--variant",
                variant,
                "--lengths",
                str(length),
                "--beam",
                str(args.beam),
                "--repetitions",
                str(args.repetitions),
                "--warmup",
                str(args.warmup),
                "--seed",
                str(args.seed),
                "--output",
                str(row_path),
            ]
            wrapped = run_command(command, args.cpu)
            print("RUN", " ".join(wrapped), flush=True)
            subprocess.run(wrapped, cwd=ROOT, check=True)
            row = read_one_csv(row_path)
            row.update(
                {
                    "side": side,
                    "executable": str(executable),
                    "executable_sha256": sha256_file(executable),
                    "invocation": str(invocation),
                }
            )
            rows.append(row)
    fieldnames = list(rows[0])
    with (root / "micro_rows.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    write_json(
        root / "run_context.json",
        {
            "schema_version": 1,
            "created_utc": utc_now(),
            "side": side,
            "cpu": args.cpu,
            "warmup": args.warmup,
            "repetitions": args.repetitions,
            "beam": args.beam,
            "seed": args.seed,
            "lengths": list(lengths),
            "variants": list(variants),
            "host": host_context(args.cpu),
            "compile": str(root / "compile.json"),
            "micro_rows": str(root / "micro_rows.csv"),
        },
    )
    print(root / "micro_rows.csv")
    return 0


def float_or_none(value: str | None) -> float | None:
    if value is None or value == "" or value.lower() == "nan":
        return None
    return float(value)


def bound_valid(row: dict[str, str]) -> bool:
    if row["variant"] == "beam":
        return True
    score = float_or_none(row.get("score"))
    upper = float_or_none(row.get("upper_bound"))
    pruned = float_or_none(row.get("pruned_upper_bound"))
    additive = float_or_none(row.get("additive_upper_bound"))
    if None in (score, upper, pruned, additive):
        return False
    if not all(math.isfinite(value) for value in (score, upper, pruned, additive)):
        return False
    tolerance = 1e-5 * max(1.0, abs(upper or 0.0), abs(additive or 0.0))
    return bool(
        score <= upper + tolerance
        and upper <= pruned + tolerance
        and upper <= additive + tolerance
        and pruned >= score - tolerance
        and additive >= score - tolerance
    )


def micro_compare(args: argparse.Namespace) -> int:
    baseline_root = resolve_path(args.baseline_root)
    candidate_root = resolve_path(args.candidate_root)
    baseline = {
        (row["variant"], int(row["length"])): row
        for row in read_micro_rows(baseline_root / "micro_rows.csv")
    }
    candidate = {
        (row["variant"], int(row["length"])): row
        for row in read_micro_rows(candidate_root / "micro_rows.csv")
    }
    keys = sorted(set(baseline) & set(candidate))
    if set(baseline) != set(candidate):
        missing = sorted(set(baseline) - set(candidate))
        unexpected = sorted(set(candidate) - set(baseline))
        raise RuntimeError(f"micro selections differ: missing={missing} unexpected={unexpected}")
    rows: list[dict[str, Any]] = []
    digest_fields = ("params_digest", "input_digest", "score_digest", "structure_digest", "bounds_digest")
    for key in keys:
        old, new = baseline[key], candidate[key]
        old_wall = float(old["wall_seconds"])
        new_wall = float(new["wall_seconds"])
        old_cpu = float(old["cpu_seconds"])
        new_cpu = float(new["cpu_seconds"])
        digest_equal = {field: old.get(field) == new.get(field) for field in digest_fields}
        rows.append(
            {
                "variant": key[0],
                "length": key[1],
                "support_mode": old.get("support_mode"),
                "beam": int(old["beam"]),
                "raw_support": int(old["raw_support"]),
                "unique_support": int(old["unique_support"]),
                "baseline_score": float_or_none(old.get("score")),
                "candidate_score": float_or_none(new.get("score")),
                "score_delta": (float_or_none(new.get("score")) or 0.0)
                - (float_or_none(old.get("score")) or 0.0),
                "baseline_wall_seconds": old_wall,
                "candidate_wall_seconds": new_wall,
                "wall_speedup": old_wall / new_wall if new_wall > 0 else None,
                "baseline_cpu_seconds": old_cpu,
                "candidate_cpu_seconds": new_cpu,
                "cpu_speedup": old_cpu / new_cpu if new_cpu > 0 else None,
                "baseline_max_rss_kb": int(old["max_rss_kb"]),
                "candidate_max_rss_kb": int(new["max_rss_kb"]),
                "rss_delta_kb": int(new["max_rss_kb"]) - int(old["max_rss_kb"]),
                "digest_equal": digest_equal,
                "input_params_equal": digest_equal["params_digest"] and digest_equal["input_digest"],
                "output_equal": digest_equal["score_digest"] and digest_equal["structure_digest"],
                "baseline_bounds_valid": bound_valid(old),
                "candidate_bounds_valid": bound_valid(new),
                "baseline_stable": old.get("output_stable") == "true",
                "candidate_stable": new.get("output_stable") == "true",
            }
        )
    output_root = resolve_path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with (output_root / "micro_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    write_json(
        output_root / "micro_comparison.json",
        {
            "schema_version": 1,
            "created_utc": utc_now(),
            "baseline_root": str(baseline_root),
            "candidate_root": str(candidate_root),
            "rows": rows,
            "output_difference_is_nonfatal": True,
            "validity_gate": "input/params, per-side stability, and bounds must hold",
        },
    )
    print(output_root / "micro_comparison.csv")
    return 0 if all(
        row["input_params_equal"]
        and row["baseline_bounds_valid"]
        and row["candidate_bounds_valid"]
        and row["baseline_stable"]
        and row["candidate_stable"]
        for row in rows
    ) else 1


def config_dataset_and_condition(config_path: Path) -> tuple[Path, list[str]]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    manifests = config.get("dataset_manifests") or []
    if not manifests:
        raise ValueError("config has no dataset manifest")
    manifest_path = resolve_path(manifests[0], config_path.parent)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    dataset = next(item for item in manifest["datasets"] if item["id"] == SERIAL_DATASET)
    condition = next(item for item in config["conditions"] if item["id"] == SERIAL_CONDITION)
    return resolve_path(dataset["input"], manifest_path.parent), [
        str(argument) for argument in condition["args"]
    ]


def timed_metrics_off(
    binary: Path,
    arguments: list[str],
    input_path: Path,
    run_dir: Path,
    timeout_seconds: float,
    cpu: int | None,
) -> dict[str, Any]:
    run_dir.mkdir(parents=True, exist_ok=True)
    prediction = run_dir / "prediction.fa"
    stderr = run_dir / "stderr.log"
    resource = run_dir / "resource.json"
    command = [str(binary), *arguments, str(input_path)]
    time_binary = Path("/usr/bin/time")
    wrapped = [
        str(time_binary),
        "-f",
        '{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kb":%M,"exit_status":%x}',
        "-o",
        str(resource),
        "--",
        *run_command(command, cpu),
    ]
    with prediction.open("wb") as output, stderr.open("wb") as error:
        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = "1"
        try:
            completed = subprocess.run(
                wrapped, cwd=ROOT, stdout=output, stderr=error,
                timeout=timeout_seconds, check=False, env=environment
            )
            status = "success" if completed.returncode == 0 else "failed"
            return_code: int | None = completed.returncode
        except subprocess.TimeoutExpired:
            status = "timeout"
            return_code = -9
    resource_value = json.loads(resource.read_text(encoding="utf-8")) if resource.is_file() else {}
    return {
        "status": status,
        "return_code": return_code,
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "command": command,
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "resource": resource_value,
        "prediction_sha256": sha256_file(prediction) if prediction.is_file() else None,
        "prediction_path": str(prediction),
    }


def prediction_diff(old_path: Path, new_path: Path) -> dict[str, Any]:
    old = old_path.read_bytes() if old_path.is_file() else b""
    new = new_path.read_bytes() if new_path.is_file() else b""
    overlap = min(len(old), len(new))
    differing = sum(left != right for left, right in zip(old[:overlap], new[:overlap]))
    return {
        "equal": old == new,
        "baseline_bytes": len(old),
        "candidate_bytes": len(new),
        "byte_diff_count": differing + abs(len(old) - len(new)),
        "byte_diff_fraction": (
            (differing + abs(len(old) - len(new))) / max(len(old), len(new))
            if max(len(old), len(new)) else 0.0
        ),
    }


def run_serial(args: argparse.Namespace) -> int:
    config_path = resolve_path(args.config)
    input_path, condition_args = config_dataset_and_condition(config_path)
    baseline = resolve_path(args.baseline_binary)
    candidate = resolve_path(args.candidate_binary)
    if sha256_file(baseline) != args.baseline_sha256:
        raise RuntimeError(f"baseline SHA256 mismatch: {baseline}")
    candidate_only_args = set(args.candidate_only_arg)
    baseline_condition_args = [
        argument for argument in condition_args if argument not in candidate_only_args
    ]
    if any(condition_args.count(argument) != 1 for argument in candidate_only_args):
        raise RuntimeError(
            "each --candidate-only-arg must occur exactly once in the candidate condition"
        )
    root = resolve_path(args.results_root)
    ensure_fresh(root)
    rows: list[dict[str, Any]] = []
    for repetition in range(1, args.repetitions + 1):
        order = (
            [("baseline", baseline), ("candidate", candidate)]
            if repetition % 2
            else [("candidate", candidate), ("baseline", baseline)]
        )
        results: dict[str, dict[str, Any]] = {}
        for order_index, (side, binary) in enumerate(order, 1):
            run_arguments = (
                condition_args if side == "candidate" else baseline_condition_args
            )
            result = timed_metrics_off(
                binary, run_arguments, input_path,
                root / side / f"rep-{repetition:02d}",
                args.timeout_seconds, args.cpu,
            )
            result.update({"side": side, "repetition": repetition, "order_index": order_index})
            write_json(root / side / f"rep-{repetition:02d}" / "result.json", result)
            results[side] = result
        old, new = results["baseline"], results["candidate"]
        old_resource, new_resource = old["resource"], new["resource"]
        old_wall, new_wall = float(old_resource.get("wall_seconds", "nan")), float(new_resource.get("wall_seconds", "nan"))
        old_cpu = float(old_resource.get("user_seconds", 0)) + float(old_resource.get("system_seconds", 0))
        new_cpu = float(new_resource.get("user_seconds", 0)) + float(new_resource.get("system_seconds", 0))
        diff = prediction_diff(Path(old["prediction_path"]), Path(new["prediction_path"]))
        rows.append(
            {
                "repetition": repetition,
                "order": "baseline-candidate" if repetition % 2 else "candidate-baseline",
                "baseline_status": old["status"],
                "candidate_status": new["status"],
                "baseline_wall_seconds": old_wall,
                "candidate_wall_seconds": new_wall,
                "wall_speedup": old_wall / new_wall if new_wall > 0 else None,
                "baseline_cpu_seconds": old_cpu,
                "candidate_cpu_seconds": new_cpu,
                "cpu_speedup": old_cpu / new_cpu if new_cpu > 0 else None,
                "baseline_max_rss_kb": old_resource.get("max_rss_kb"),
                "candidate_max_rss_kb": new_resource.get("max_rss_kb"),
                "prediction_bytes_equal": diff["equal"],
                "prediction_byte_diff_count": diff["byte_diff_count"],
                "prediction_byte_diff_fraction": diff["byte_diff_fraction"],
            }
        )
    with (root / "paired.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_json(
        root / "paired.json",
        {
            "schema_version": 1,
            "created_utc": utc_now(),
            "config": str(config_path),
            "dataset": SERIAL_DATASET,
            "condition": SERIAL_CONDITION,
            "metrics_enabled": False,
            "baseline_sha256": args.baseline_sha256,
            "candidate_only_args": sorted(candidate_only_args),
            "cpu": args.cpu,
            "repetitions": args.repetitions,
            "rows": rows,
        },
    )
    print(root / "paired.csv")
    return 0 if all(
        row["baseline_status"] == "success" and row["candidate_status"] == "success"
        for row in rows
    ) else 1


def result_paths(root: Path) -> dict[tuple[str, str, int], Path]:
    result: dict[tuple[str, str, int], Path] = {}
    for path in sorted((root / "runs").glob("*/*/rep-*/result.json")):
        data = json.loads(path.read_text(encoding="utf-8"))
        key = (str(data["condition"]), str(data["dataset"]), int(data["repetition"]))
        if key in result:
            raise RuntimeError(f"duplicate result key: {key}")
        result[key] = path
    return result


def load_events(result: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    path = Path(result.get("artifacts", {}).get("metrics", ""))
    if not path.is_file():
        raise RuntimeError(f"missing metrics trace: {path}")
    events, validation = runner_load_metrics(path)
    if validation.get("parse_errors"):
        raise RuntimeError(
            f"malformed metrics trace {path}: {validation['parse_errors']}"
        )
    if not validation.get("run_summary"):
        raise RuntimeError(f"metrics trace has no run_summary: {path}")
    return events, validation


def trace_metrics(result: dict[str, Any]) -> dict[str, Any]:
    events, validation = load_events(result)
    summary = validation.get("run_summary", {}) or {}
    dd_summaries = [event for event in events if event.get("event") == "dd_summary"]
    dd_iterations = [event for event in events if event.get("event") == "dd_iteration"]
    violations = [float(event["violated"]) for event in dd_summaries if event.get("violated") is not None]
    last_iteration_by_merge: dict[str, dict[str, Any]] = {}
    for event in dd_iterations:
        last_iteration_by_merge[str(event["merge_id"])] = event
    certified = [
        float(event["certified_track_violated"])
        for event in last_iteration_by_merge.values()
        if event.get("certified_track_violated") is not None
    ]
    gaps = [float(event["gap"]) for event in dd_summaries if event.get("gap") is not None]
    resource = result.get("resource", {}) or {}
    quality = result.get("quality", {}) or {}
    return {
        "status": result.get("status"),
        "metrics_valid": bool(validation.get("valid")),
        "parse_error_count": len(validation.get("parse_errors", []) or []),
        "invariant_error_count": len(validation.get("invariant_errors", []) or []),
        "mcc": quality.get("mcc"),
        "sci": quality.get("sci"),
        "sps": quality.get("sps"),
        "dd_iterations": summary.get("dd_iterations"),
        "dd_calls": summary.get("dd_calls"),
        "dd_seconds": summary.get("dd_seconds"),
        "wall_seconds": resource.get("wall_seconds"),
        "cpu_seconds": (
            float(resource.get("user_seconds", 0))
            + float(resource.get("system_seconds", 0))
        ),
        "max_rss_kb": resource.get("max_rss_kb"),
        "dd_violations_total": sum(violations) if violations else 0.0,
        "dd_violations_mean_merge": statistics.fmean(violations) if violations else 0.0,
        "dd_violations_max_merge": max(violations) if violations else 0.0,
        "certified_violations_total": sum(certified) if certified else 0.0,
        "certified_violations_mean_merge": statistics.fmean(certified) if certified else 0.0,
        "certified_violations_max_merge": max(certified) if certified else 0.0,
        "dd_final_gap_mean": statistics.fmean(gaps) if gaps else 0.0,
        "dd_final_gap_max": max(gaps) if gaps else 0.0,
    }


def numeric_delta(old: Any, new: Any) -> float | None:
    if old is None or new is None:
        return None
    return float(new) - float(old)


def manifest_selection_digest(manifest: dict[str, Any]) -> str:
    config = manifest.get("config") or {}
    selection = {
        "dataset_manifests": config.get("dataset_manifests"),
        "dataset_filter": config.get("dataset_filter"),
        "conditions": config.get("conditions"),
        "repetitions": config.get("repetitions"),
        "random_seed": config.get("random_seed"),
        "environment": config.get("environment"),
        "sci": config.get("sci"),
    }
    payload = json.dumps(selection, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def load_result_manifest(root: Path) -> dict[str, Any]:
    path = root / "manifest.json"
    if not path.is_file():
        raise RuntimeError(f"missing result manifest: {path}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    for field in ("config_sha256", "binary_sha256"):
        if not manifest.get(field):
            raise RuntimeError(f"manifest missing {field}: {path}")
    manifest["selection_digest"] = manifest_selection_digest(manifest)
    return manifest


def normalized_command(result: dict[str, Any]) -> list[str]:
    command = list(result.get("command") or [])
    if not command:
        return []
    normalized = ["<binary>"]
    index = 1
    while index < len(command):
        argument = command[index]
        if argument == "--metrics-jsonl":
            index += 2
            continue
        if index == len(command) - 1:
            normalized.append("<input>")
        else:
            normalized.append(argument)
        index += 1
    return normalized


def compare_full(args: argparse.Namespace) -> int:
    baseline_root = resolve_path(args.baseline_results)
    candidate_root = resolve_path(args.candidate_results)
    baseline_manifest = load_result_manifest(baseline_root)
    candidate_manifest = load_result_manifest(candidate_root)
    baseline_paths = result_paths(baseline_root)
    candidate_paths = result_paths(candidate_root)
    if not baseline_paths or not candidate_paths:
        raise RuntimeError("full comparison requires non-empty baseline and candidate result sets")
    if len(candidate_paths) != args.expected_runs:
        raise RuntimeError(
            f"candidate result count {len(candidate_paths)} != expected {args.expected_runs}"
        )
    if args.allow_baseline_superset:
        if not set(candidate_paths).issubset(baseline_paths):
            raise RuntimeError("candidate result keys are not a subset of baseline keys")
    elif len(baseline_paths) != args.expected_runs:
        raise RuntimeError(
            f"baseline result count {len(baseline_paths)} != expected {args.expected_runs}"
        )
    common = sorted(set(baseline_paths) & set(candidate_paths))
    missing = sorted(set(baseline_paths) - set(candidate_paths))
    unexpected = sorted(set(candidate_paths) - set(baseline_paths))
    if missing or unexpected:
        print(f"run selection differs: missing={len(missing)} unexpected={len(unexpected)}", file=sys.stderr)
    rows: list[dict[str, Any]] = []
    for key in common:
        old = json.loads(baseline_paths[key].read_text(encoding="utf-8"))
        new = json.loads(candidate_paths[key].read_text(encoding="utf-8"))
        old_trace, new_trace = trace_metrics(old), trace_metrics(new)
        old_prediction = Path(old.get("artifacts", {}).get("prediction", ""))
        new_prediction = Path(new.get("artifacts", {}).get("prediction", ""))
        diff = prediction_diff(old_prediction, new_prediction)
        row: dict[str, Any] = {
            "condition": key[0],
            "dataset": key[1],
            "repetition": key[2],
            "baseline_status": old_trace["status"],
            "candidate_status": new_trace["status"],
            "input_sha256_equal": bool(old.get("input_sha256"))
            and old.get("input_sha256") == new.get("input_sha256"),
            "condition_args_equal": normalized_command(old) == normalized_command(new),
            **{f"baseline_{field}": value for field, value in old_trace.items() if field not in ("status",)},
            **{f"candidate_{field}": value for field, value in new_trace.items() if field not in ("status",)},
            "wall_speedup": (
                float(old_trace["wall_seconds"]) / float(new_trace["wall_seconds"])
                if old_trace["wall_seconds"] and new_trace["wall_seconds"] else None
            ),
            "cpu_speedup": (
                float(old_trace["cpu_seconds"]) / float(new_trace["cpu_seconds"])
                if old_trace["cpu_seconds"] and new_trace["cpu_seconds"] else None
            ),
            "mcc_delta": numeric_delta(old_trace["mcc"], new_trace["mcc"]),
            "sci_delta": numeric_delta(old_trace["sci"], new_trace["sci"]),
            "sps_delta": numeric_delta(old_trace["sps"], new_trace["sps"]),
            "dd_iterations_delta": numeric_delta(old_trace["dd_iterations"], new_trace["dd_iterations"]),
            "dd_violations_total_delta": numeric_delta(old_trace["dd_violations_total"], new_trace["dd_violations_total"]),
            "dd_violations_mean_merge_delta": numeric_delta(old_trace["dd_violations_mean_merge"], new_trace["dd_violations_mean_merge"]),
            "dd_final_gap_mean_delta": numeric_delta(old_trace["dd_final_gap_mean"], new_trace["dd_final_gap_mean"]),
            "max_rss_delta_kb": numeric_delta(old_trace["max_rss_kb"], new_trace["max_rss_kb"]),
            "prediction_bytes_equal": diff["equal"],
            "prediction_baseline_bytes": diff["baseline_bytes"],
            "prediction_candidate_bytes": diff["candidate_bytes"],
            "prediction_byte_diff_count": diff["byte_diff_count"],
            "prediction_byte_diff_fraction": diff["byte_diff_fraction"],
            "prediction_sha256_equal": old.get("artifacts", {}).get("prediction")
            and new.get("artifacts", {}).get("prediction")
            and sha256_file(old_prediction) == sha256_file(new_prediction),
        }
        rows.append(row)
    output_root = resolve_path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else ["condition", "dataset", "repetition"]
    with (output_root / "paired.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(row["condition"], []).append(row)
    condition_rows: list[dict[str, Any]] = []
    numeric_fields = (
        "wall_speedup", "cpu_speedup", "mcc_delta", "sci_delta", "sps_delta",
        "dd_iterations_delta", "dd_violations_total_delta", "dd_violations_mean_merge_delta",
        "dd_final_gap_mean_delta", "max_rss_delta_kb", "prediction_byte_diff_fraction",
        "baseline_wall_seconds", "candidate_wall_seconds", "baseline_cpu_seconds",
        "candidate_cpu_seconds", "baseline_max_rss_kb", "candidate_max_rss_kb",
        "baseline_mcc", "candidate_mcc", "baseline_sci", "candidate_sci",
        "baseline_sps", "candidate_sps", "baseline_dd_iterations", "candidate_dd_iterations",
        "baseline_dd_violations_total", "candidate_dd_violations_total",
    )
    for condition, condition_values in sorted(grouped.items()):
        row: dict[str, Any] = {"condition": condition, "n": len(condition_values)}
        for field in numeric_fields:
            values = [float(value[field]) for value in condition_values if value.get(field) is not None]
            row[f"mean_{field}"] = statistics.fmean(values) if values else None
        for field in (
            "baseline_metrics_valid", "candidate_metrics_valid",
            "prediction_bytes_equal", "baseline_status", "candidate_status",
        ):
            values = [value[field] for value in condition_values]
            if field.endswith("metrics_valid") or field == "prediction_bytes_equal":
                row[f"count_{field}_true"] = sum(value is True for value in values)
            else:
                row[f"count_{field}_success"] = sum(value == "success" for value in values)
        condition_rows.append(row)
    with (output_root / "per_condition.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(condition_rows[0]) if condition_rows else ["condition", "n"])
        writer.writeheader()
        writer.writerows(condition_rows)

    invalid_rows = [
        row for row in rows
        if not row["baseline_metrics_valid"] or not row["candidate_metrics_valid"]
        or row["baseline_status"] != "success" or row["candidate_status"] != "success"
        or not row["input_sha256_equal"]
        or not row["condition_args_equal"]
    ]
    summary = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "baseline_results": str(baseline_root),
        "candidate_results": str(candidate_root),
        "expected_runs": args.expected_runs,
        "expected_common_runs": len(common),
        "baseline_manifest_config_sha256": baseline_manifest["config_sha256"],
        "candidate_manifest_config_sha256": candidate_manifest["config_sha256"],
        "baseline_manifest_binary_sha256": baseline_manifest["binary_sha256"],
        "candidate_manifest_binary_sha256": candidate_manifest["binary_sha256"],
        "baseline_selection_digest": baseline_manifest["selection_digest"],
        "candidate_selection_digest": candidate_manifest["selection_digest"],
        "selection_digest_equal": baseline_manifest["selection_digest"] == candidate_manifest["selection_digest"],
        "missing_runs": [list(value) for value in missing],
        "unexpected_runs": [list(value) for value in unexpected],
        "rows": len(rows),
        "invalid_or_failed_rows": len(invalid_rows),
        "input_sha256_equal_count": sum(row["input_sha256_equal"] is True for row in rows),
        "condition_args_equal_count": sum(row["condition_args_equal"] is True for row in rows),
        "output_difference_is_nonfatal": True,
        "validity_gate": "success status and metrics_validation.valid for every paired run",
        "quality_fields": ["mcc", "sci", "sps"],
        "trace_fields": ["dd_iterations", "dd_violations_total", "dd_final_gap_mean"],
        "per_condition": condition_rows,
    }
    write_json(output_root / "comparison.json", summary)
    print(output_root / "per_condition.csv")
    return 0 if (
        len(common) == args.expected_runs
        and not missing
        and not unexpected
        and not invalid_rows
        and summary["selection_digest_equal"]
        and all(row["input_sha256_equal"] and row["condition_args_equal"] for row in rows)
    ) else 1


def add_micro_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--lengths", type=int, nargs="+", default=list(DEFAULT_LENGTHS))
    parser.add_argument("--variants", nargs="+", default=list(DEFAULT_VARIANTS))
    parser.add_argument("--beam", type=int, default=100)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--repetitions", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cpu", type=int, default=DEFAULT_CPU)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    compile_parser = subparsers.add_parser("micro-compile")
    compile_parser.add_argument("--side", choices=("baseline", "candidate"), required=True)
    compile_parser.add_argument("--root", type=Path, required=True)
    compile_parser.add_argument("--source-root", type=Path)
    compile_parser.set_defaults(handler=compile_micro)

    run_parser = subparsers.add_parser("micro-run")
    add_micro_options(run_parser)
    run_parser.set_defaults(handler=run_micro)

    compare_micro_parser = subparsers.add_parser("micro-compare")
    compare_micro_parser.add_argument("--baseline-root", type=Path, required=True)
    compare_micro_parser.add_argument("--candidate-root", type=Path, required=True)
    compare_micro_parser.add_argument("--output-root", type=Path, required=True)
    compare_micro_parser.set_defaults(handler=micro_compare)

    serial_parser = subparsers.add_parser("serial", help="metrics-off paired AB")
    serial_parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    serial_parser.add_argument("--baseline-binary", type=Path, default=DEFAULT_BASELINE_BINARY)
    serial_parser.add_argument("--candidate-binary", type=Path, default=ROOT / "build-arm64/src/dafs")
    serial_parser.add_argument("--baseline-sha256", default=BASELINE_SHA256)
    serial_parser.add_argument(
        "--candidate-only-arg",
        action="append",
        default=[],
        help="remove this exact option from the old-baseline invocation only",
    )
    serial_parser.add_argument("--results-root", type=Path, required=True)
    serial_parser.add_argument("--repetitions", type=int, default=5)
    serial_parser.add_argument("--timeout-seconds", type=float, default=3600.0)
    serial_parser.add_argument("--cpu", type=int, default=DEFAULT_CPU)
    serial_parser.set_defaults(handler=run_serial)

    full_parser = subparsers.add_parser("compare", help="approximate full-result comparison")
    full_parser.add_argument("--baseline-results", type=Path, default=DEFAULT_BASELINE_RESULTS)
    full_parser.add_argument("--candidate-results", type=Path, default=DEFAULT_CANDIDATE_RESULTS)
    full_parser.add_argument("--output-root", type=Path, required=True)
    full_parser.add_argument("--expected-runs", type=int, default=680)
    full_parser.add_argument("--allow-baseline-superset", action="store_true")
    full_parser.set_defaults(handler=compare_full)

    args = parser.parse_args()
    try:
        return int(args.handler(args))
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
