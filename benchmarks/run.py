#!/usr/bin/env python3
"""Reproducible, resumable benchmark runner for DAFS."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import random
import re
import shlex
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from score import parse_alignment, score


SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True,
                                    allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(root: Path, *arguments: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *arguments], cwd=root, text=True,
            stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def load_resource(path: Path) -> dict[str, Any]:
    if not path.exists() or not path.read_text(encoding="utf-8").strip():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        return {"parse_error": str(error), "raw": path.read_text(encoding="utf-8")}


def load_metrics(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    events: list[dict[str, Any]] = []
    errors: list[str] = []
    if path.exists():
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            try:
                events.append(json.loads(line))
            except json.JSONDecodeError as error:
                errors.append(f"line {number}: {error}")

    by_merge: dict[int, list[dict[str, Any]]] = {}
    for event in events:
        if event.get("event") == "dd_iteration":
            by_merge.setdefault(int(event["merge_id"]), []).append(event)
    invariant_errors: list[str] = []
    for merge_id, iterations in by_merge.items():
        previous_ub = float("inf")
        previous_lb = float("-inf")
        for event in iterations:
            ub, lb = event.get("best_ub"), event.get("lb")
            if ub is not None and lb is not None and lb > ub + 1e-5 * max(1.0, abs(ub)):
                invariant_errors.append(f"merge {merge_id}: LB {lb} > UB {ub}")
            if ub is not None and ub > previous_ub + 1e-5 * max(1.0, abs(previous_ub)):
                invariant_errors.append(f"merge {merge_id}: BestUB increased")
            if lb is not None and lb < previous_lb - 1e-5 * max(1.0, abs(previous_lb)):
                invariant_errors.append(f"merge {merge_id}: LB decreased")
            if ub is not None:
                previous_ub = ub
            if lb is not None:
                previous_lb = lb
    summary = next((event for event in reversed(events)
                    if event.get("event") == "run_summary"), {})
    validation = {
        "event_count": len(events),
        "parse_errors": errors,
        "invariant_errors": invariant_errors,
        "valid": not errors and not invariant_errors and bool(summary),
        "run_summary": summary,
    }
    return events, validation


def resolve_path(value: str, base: Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def validate_id(value: str, kind: str) -> str:
    if not SAFE_ID.fullmatch(value):
        raise ValueError(f"{kind} id must match {SAFE_ID.pattern}: {value!r}")
    return value


def command_digest(command: list[str], input_sha256: str) -> str:
    payload = json.dumps({"command": command, "input_sha256": input_sha256},
                         sort_keys=True).encode()
    return hashlib.sha256(payload).hexdigest()


def execute(command: list[str], stdout_path: Path, stderr_path: Path,
            timeout_seconds: float, environment: dict[str, str]) -> tuple[int | None, bool]:
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                                   env=environment, start_new_session=True)
        try:
            return process.wait(timeout=timeout_seconds), False
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            return None, True


def flatten_summary(result: dict[str, Any]) -> dict[str, Any]:
    resource = result.get("resource", {})
    quality = result.get("quality", {}) or {}
    internal = result.get("metrics_validation", {}).get("run_summary", {})
    return {
        "condition": result["condition"], "dataset": result["dataset"],
        "repetition": result["repetition"], "status": result["status"],
        "return_code": result.get("return_code"),
        "wall_seconds": resource.get("wall_seconds"),
        "user_seconds": resource.get("user_seconds"),
        "system_seconds": resource.get("system_seconds"),
        "max_rss_kb": resource.get("max_rss_kb"),
        "sps": quality.get("sps"), "sensitivity": quality.get("sensitivity"),
        "ppv": quality.get("ppv"), "mcc": quality.get("mcc"),
        "cbp_f1": quality.get("cbp_f1"),
        "quality_complete": result.get("quality_complete"),
        "internal_seconds": internal.get("seconds"),
        "dd_iterations": internal.get("dd_iterations"),
        "cbp_peak": internal.get("cbp_peak"),
        "metrics_valid": result.get("metrics_validation", {}).get("valid"),
    }


def write_summaries(output_dir: Path) -> None:
    results = []
    for path in sorted((output_dir / "runs").glob("*/*/rep-*/result.json")):
        results.append(json.loads(path.read_text(encoding="utf-8")))
    with (output_dir / "summary.jsonl").open("w", encoding="utf-8") as stream:
        for result in results:
            stream.write(json.dumps(flatten_summary(result), sort_keys=True,
                                    allow_nan=False) + "\n")
    rows = [flatten_summary(result) for result in results]
    fieldnames = list(rows[0]) if rows else list(flatten_summary({
        "condition": "", "dataset": "", "repetition": 0, "status": ""
    }))
    with (output_dir / "summary.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--binary", type=Path, help="override config binary")
    parser.add_argument("--output-dir", type=Path, help="override config output_dir")
    parser.add_argument("--rerun", action="store_true", help="rerun completed cases")
    parser.add_argument("--rerun-failed", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    config_path = args.config.resolve()
    config_base = config_path.parent
    config = json.loads(config_path.read_text(encoding="utf-8"))
    repository = Path(__file__).resolve().parents[1]
    binary = (args.binary.resolve() if args.binary else
              resolve_path(config["binary"], config_base))
    output_dir = (args.output_dir.resolve() if args.output_dir else
                  resolve_path(config.get("output_dir", "results"), config_base))
    if not binary.is_file():
        raise FileNotFoundError(f"DAFS binary not found: {binary}")
    time_binary = Path(config.get("time_binary", "/usr/bin/time"))
    if not time_binary.is_file():
        raise FileNotFoundError(f"GNU time not found: {time_binary}")

    repetitions = int(config.get("repetitions", 1))
    timeout_seconds = float(config.get("timeout_seconds", 3600))
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    environment.update({str(key): str(value)
                        for key, value in config.get("environment", {}).items()})
    datasets = list(config.get("datasets", []))
    for manifest_name in config.get("dataset_manifests", []):
        manifest_path = resolve_path(manifest_name, config_base)
        manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest_base = manifest_path.parent
        for dataset in manifest_data["datasets"]:
            resolved = dict(dataset)
            resolved["input"] = str(resolve_path(dataset["input"], manifest_base))
            if dataset.get("reference"):
                resolved["reference"] = str(resolve_path(dataset["reference"], manifest_base))
            datasets.append(resolved)
    if not datasets:
        raise ValueError("config must define datasets or dataset_manifests")
    conditions = config["conditions"]
    jobs = [(dataset, condition, repetition)
            for dataset in datasets for condition in conditions
            for repetition in range(1, repetitions + 1)]
    random.Random(int(config.get("random_seed", 42))).shuffle(jobs)

    if not args.dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "runs").mkdir(exist_ok=True)
        manifest = {
            "schema_version": 1,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "config": config,
            "config_path": str(config_path),
            "config_sha256": sha256_file(config_path),
            "binary": str(binary),
            "binary_sha256": sha256_file(binary),
            "git_commit": git_value(repository, "rev-parse", "HEAD"),
            "git_status": git_value(repository, "status", "--short"),
            "platform": platform.platform(),
            "uname": list(platform.uname()),
            "python": sys.version,
            "environment": {key: environment[key] for key in sorted(config.get("environment", {}))},
            "omp_num_threads": environment.get("OMP_NUM_THREADS"),
        }
        atomic_json(output_dir / "manifest.json", manifest)

    for dataset, condition, repetition in jobs:
        dataset_id = validate_id(str(dataset["id"]), "dataset")
        condition_id = validate_id(str(condition["id"]), "condition")
        input_path = resolve_path(dataset["input"], config_base)
        reference_path = (resolve_path(dataset["reference"], config_base)
                          if dataset.get("reference") else None)
        if not input_path.is_file():
            raise FileNotFoundError(f"input not found: {input_path}")
        if reference_path is not None and not reference_path.is_file():
            raise FileNotFoundError(f"reference not found: {reference_path}")
        run_dir = output_dir / "runs" / condition_id / dataset_id / f"rep-{repetition:02d}"
        metrics_path = run_dir / "metrics.jsonl"
        prediction_path = run_dir / "prediction.fa"
        stderr_path = run_dir / "stderr.log"
        resource_path = run_dir / "resource.json"
        result_path = run_dir / "result.json"
        dafs_command = [str(binary), *map(str, condition.get("args", [])),
                        "--metrics-jsonl", str(metrics_path), str(input_path)]
        digest = command_digest(dafs_command, sha256_file(input_path))
        if result_path.exists() and not args.rerun:
            previous = json.loads(result_path.read_text(encoding="utf-8"))
            failed = previous.get("status") != "success"
            if previous.get("command_sha256") == digest and not (args.rerun_failed and failed):
                print(f"SKIP {condition_id}/{dataset_id}/rep-{repetition:02d}")
                continue
        print("RUN ", shlex.join(dafs_command))
        if args.dry_run:
            continue
        run_dir.mkdir(parents=True, exist_ok=True)
        time_format = ("{\"wall_seconds\":%e,\"user_seconds\":%U,"
                       "\"system_seconds\":%S,\"max_rss_kb\":%M,"
                       "\"exit_status\":%x}")
        wrapped = [str(time_binary), "-f", time_format, "-o", str(resource_path),
                   "--", *dafs_command]
        started = datetime.now(timezone.utc)
        monotonic_start = time.monotonic()
        return_code, timed_out = execute(wrapped, prediction_path, stderr_path,
                                         timeout_seconds, environment)
        harness_seconds = time.monotonic() - monotonic_start
        _, metrics_validation = load_metrics(metrics_path)
        quality = None
        score_error = None
        if return_code == 0 and reference_path is not None:
            try:
                quality = score(parse_alignment(prediction_path),
                                parse_alignment(reference_path))
            except Exception as error:  # retain failed scoring as benchmark evidence
                score_error = f"{type(error).__name__}: {error}"
        result = {
            "schema_version": 1,
            "condition": condition_id,
            "dataset": dataset_id,
            "repetition": repetition,
            "status": "timeout" if timed_out else ("success" if return_code == 0 else "failed"),
            "return_code": return_code,
            "timed_out": timed_out,
            "started_utc": started.isoformat(),
            "harness_wall_seconds": harness_seconds,
            "command": dafs_command,
            "command_sha256": digest,
            "input": str(input_path),
            "input_sha256": sha256_file(input_path),
            "reference": str(reference_path) if reference_path else None,
            "resource": load_resource(resource_path),
            "metrics_validation": metrics_validation,
            "quality": quality,
            "quality_complete": (None if reference_path is None else
                                 bool(quality is not None and
                                      not quality["missing_predicted_sequences"] and
                                      not quality["extra_predicted_sequences"])),
            "score_error": score_error,
            "artifacts": {
                "prediction": str(prediction_path), "stderr": str(stderr_path),
                "metrics": str(metrics_path), "resource": str(resource_path),
            },
        }
        atomic_json(result_path, result)
        print(f"{result['status'].upper()} {condition_id}/{dataset_id}/rep-{repetition:02d}")

    if not args.dry_run:
        write_summaries(output_dir)
        print(f"Summary: {output_dir / 'summary.csv'}")
        failures = []
        for path in sorted((output_dir / "runs").glob("*/*/rep-*/result.json")):
            result = json.loads(path.read_text(encoding="utf-8"))
            if result.get("status") != "success":
                failures.append(f"{path}: status={result.get('status')}")
            elif not result.get("metrics_validation", {}).get("valid"):
                failures.append(f"{path}: invalid internal metrics")
            elif result.get("reference") and result.get("score_error"):
                failures.append(f"{path}: {result['score_error']}")
            elif result.get("reference") and not result.get("quality_complete"):
                failures.append(f"{path}: prediction/reference sequence sets differ")
        if failures and not config.get("allow_failures", False):
            for failure in failures:
                print(f"ERROR {failure}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
