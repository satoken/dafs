#!/usr/bin/env python3
"""Run LinearTurboFold on the benchmark manifests used by DAFS.

The official LinearTurboFold binary writes an aligned FASTA file and one
dot-bracket ``.db`` file per sequence.  This runner keeps those artifacts,
projects the individual structures onto the alignment, and evaluates them
with the common benchmark scorer.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import platform
import random
import shlex
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from run import (atomic_json, command_digest, execute,
                 git_value, load_resource, resolve_path, resolve_sci_config,
                 resolve_time_binary, sha256_file, validate_id, write_summaries)
from score import (calculate_sci, match_prediction_names, parse_alignment,
                   parse_linearturbofold, score)


def load_datasets(config: dict[str, Any], config_base: Path) -> list[dict[str, Any]]:
    datasets = [dict(dataset) for dataset in config.get("datasets", [])]
    for manifest_name in config.get("dataset_manifests", []):
        manifest_path = resolve_path(manifest_name, config_base)
        manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest_base = manifest_path.parent
        for dataset in manifest_data["datasets"]:
            resolved = dict(dataset)
            resolved["input"] = str(resolve_path(dataset["input"], manifest_base))
            if dataset.get("reference"):
                resolved["reference"] = str(
                    resolve_path(dataset["reference"], manifest_base))
            datasets.append(resolved)

    dataset_filter = config.get("dataset_filter", {})
    if dataset_filter.get("collections"):
        collections = set(map(str, dataset_filter["collections"]))
        datasets = [dataset for dataset in datasets
                    if str(dataset.get("collection")) in collections]
    if dataset_filter.get("ids"):
        dataset_ids = set(map(str, dataset_filter["ids"]))
        datasets = [dataset for dataset in datasets
                    if str(dataset.get("id")) in dataset_ids]
    if not datasets:
        raise ValueError("config and dataset_filter selected no datasets")
    return datasets


def option_values(condition: dict[str, Any]) -> list[str]:
    """Convert a readable condition object to the binary's positional API."""
    if "args" in condition:
        values = [str(value) for value in condition["args"]]
        if len(values) != 9:
            raise ValueError(
                "LinearTurboFold condition.args must contain exactly 9 values: "
                "hmm_beam cky_beam iterations save_bpps save_pfs verbose "
                "min_helix_length pk_iterations threshold")
        return values

    options = {
        "hmm_beam": 100,
        "cky_beam": 100,
        "iterations": 3,
        "save_bpps": 0,
        "save_pfs": 0,
        "verbose": 0,
        "min_helix_length": 3,
        "pk_iterations": 1,
        "threshold": 0.3,
    }
    options.update(condition.get("options", {}))
    boolean_names = {"save_bpps", "save_pfs", "verbose"}
    values = []
    for name in ("hmm_beam", "cky_beam", "iterations", "save_bpps", "save_pfs",
                 "verbose", "min_helix_length", "pk_iterations", "threshold"):
        value = options[name]
        if name in boolean_names:
            value = int(bool(value)) if isinstance(value, bool) else int(value)
        values.append(str(value))
    return values


def clear_previous_outputs(run_dir: Path) -> None:
    """Remove only known outputs from this one resumable case directory."""
    for pattern in ("*.db", "*.ct", "*.bpp", "*.pfs", "output.aln"):
        for path in run_dir.glob(pattern):
            path.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--binary", type=Path,
                        help="override config binary (the compiled bin/linearturbofold)")
    parser.add_argument("--output-dir", type=Path, help="override config output_dir")
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--rerun-failed", action="store_true")
    parser.add_argument(
        "--score-only", action="store_true",
        help="score existing output.aln/.db artifacts without running the binary")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--sci-bin-dir", type=Path,
                        help="override SCI ViennaRNA bin_dir from config")
    parser.add_argument("--sci-version",
                        help="override the expected ViennaRNA version")
    args = parser.parse_args()
    if args.jobs <= 0:
        parser.error("--jobs must be positive")

    config_path = args.config.resolve()
    config_base = config_path.parent
    config = json.loads(config_path.read_text(encoding="utf-8"))
    sci_value = config.get("sci")
    if args.sci_bin_dir or args.sci_version:
        sci_value = dict(sci_value or {})
        if args.sci_bin_dir:
            sci_value["bin_dir"] = str(args.sci_bin_dir)
        if args.sci_version:
            sci_value["version"] = args.sci_version
    sci_config = resolve_sci_config(sci_value, config_base)

    repository = Path(__file__).resolve().parents[1]
    binary = (args.binary.resolve() if args.binary else
              resolve_path(config["binary"], config_base))
    output_dir = (args.output_dir.resolve() if args.output_dir else
                  resolve_path(config.get("output_dir", "results"), config_base))
    if (not args.score_only and
            (not binary.is_file() or not os.access(binary, os.X_OK))):
        raise FileNotFoundError(
            f"LinearTurboFold binary not found or not executable: {binary}")
    time_binary = (None if args.score_only else
                   resolve_time_binary(
                       config.get("time_binary", os.environ.get("TIME_BINARY", "time"))))

    sci_provenance = (None if sci_config is None else {
        "rnaalifold": str(sci_config.rnaalifold),
        "rnafold": str(sci_config.rnafold),
        "expected_version": sci_config.expected_version,
        "timeout_seconds": sci_config.timeout_seconds,
        "rnaalifold_sha256": sha256_file(Path(sci_config.rnaalifold)),
        "rnafold_sha256": sha256_file(Path(sci_config.rnafold)),
    })
    source_dir = config.get("source_dir")
    source_path = (resolve_path(source_dir, config_base)
                   if source_dir else binary.parent.parent)
    if not args.score_only and not source_path.is_dir():
        raise FileNotFoundError(f"LinearTurboFold source directory not found: {source_path}")
    datasets = load_datasets(config, config_base)
    conditions = config.get("conditions") or [{
        "id": "linearturbofold",
        "ribosum": "not_applicable",
        "options": {},
    }]
    repetitions = int(config.get("repetitions", 1))
    timeout_seconds = float(config.get("timeout_seconds", 3600))
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    environment.update({str(key): str(value)
                        for key, value in config.get("environment", {}).items()})
    jobs = [(dataset, condition, repetition)
            for dataset in datasets for condition in conditions
            for repetition in range(1, repetitions + 1)]
    random.Random(int(config.get("random_seed", 42))).shuffle(jobs)

    if not args.dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "runs").mkdir(exist_ok=True)
        manifest = {
            "schema_version": 1,
            "method": "LinearTurboFold",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "config": config,
            "config_path": str(config_path),
            "config_sha256": sha256_file(config_path),
            "binary": str(binary),
            "binary_sha256": sha256_file(binary),
            "source_dir": str(source_path),
            "source_git_commit": git_value(source_path, "rev-parse", "HEAD"),
            "git_commit": git_value(repository, "rev-parse", "HEAD"),
            "git_status": git_value(repository, "status", "--short"),
            "platform": platform.platform(),
            "uname": list(platform.uname()),
            "python": sys.version,
            "environment": {key: environment[key]
                            for key in sorted(config.get("environment", {}))},
            "omp_num_threads": environment.get("OMP_NUM_THREADS"),
            "sci": sci_provenance,
        }
        atomic_json(output_dir / "manifest.json", manifest)

    def run_job(job: tuple[dict[str, Any], dict[str, Any], int]) -> None:
        dataset, condition, repetition = job
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
        alignment_path = run_dir / "output.aln"
        stdout_path = run_dir / "stdout.log"
        stderr_path = run_dir / "stderr.log"
        resource_path = run_dir / "resource.json"
        result_path = run_dir / "result.json"
        command = [str(binary), str(input_path), str(run_dir), *option_values(condition)]
        digest = command_digest(command, sha256_file(input_path), {
            "runner": "run_linearturbofold.py",
            "sci": sci_provenance,
        })
        if args.score_only:
            if not result_path.is_file():
                raise FileNotFoundError(
                    f"existing result.json not found for score-only run: {result_path}")
            if not alignment_path.is_file():
                raise FileNotFoundError(
                    f"LinearTurboFold alignment not found: {alignment_path}")
            previous = json.loads(result_path.read_text(encoding="utf-8"))
            quality = None
            score_error = None
            reference = None
            try:
                reference = (parse_alignment(reference_path)
                             if reference_path else None)
                prediction = parse_linearturbofold(alignment_path, run_dir)
                if reference is not None:
                    prediction = match_prediction_names(prediction, reference)
                    quality = score(prediction, reference)
                if sci_config is not None:
                    sci_quality = calculate_sci(prediction, sci_config)
                    quality = {**(quality or {}), **sci_quality}
            except Exception as error:  # preserve failed scoring as evidence
                score_error = f"{type(error).__name__}: {error}"
            previous.update({
                "quality": quality,
                "quality_complete": (None if reference is None else
                                     bool(quality is not None and
                                          not quality["missing_predicted_sequences"] and
                                          not quality["extra_predicted_sequences"])),
                "score_error": score_error,
                "score_only_utc": datetime.now(timezone.utc).isoformat(),
            })
            atomic_json(result_path, previous)
            print(f"RESCORED {condition_id}/{dataset_id}/rep-{repetition:02d}")
            return

        if result_path.exists() and not args.rerun:
            previous = json.loads(result_path.read_text(encoding="utf-8"))
            failed = previous.get("status") != "success"
            if (previous.get("command_sha256") == digest and
                    alignment_path.is_file() and
                    not (args.rerun_failed and failed)):
                print(f"SKIP {condition_id}/{dataset_id}/rep-{repetition:02d}")
                return

        print(f"RUN  {shlex.join(command)}", flush=True)
        if args.dry_run:
            return
        run_dir.mkdir(parents=True, exist_ok=True)
        clear_previous_outputs(run_dir)
        time_format = ("{\"wall_seconds\":%e,\"user_seconds\":%U,"
                       "\"system_seconds\":%S,\"max_rss_kb\":%M,"
                       "\"exit_status\":%x}")
        wrapped = [str(time_binary), "-f", time_format, "-o", str(resource_path),
                   "--", *command]
        started = datetime.now(timezone.utc)
        monotonic_start = time.monotonic()
        return_code, timed_out = execute(wrapped, stdout_path, stderr_path,
                                         timeout_seconds, environment,
                                         cwd=source_path)
        harness_seconds = time.monotonic() - monotonic_start
        quality = None
        score_error = None
        if return_code == 0:
            try:
                if not alignment_path.is_file():
                    raise FileNotFoundError(
                        f"LinearTurboFold did not produce {alignment_path}")
                prediction = parse_linearturbofold(alignment_path, run_dir)
                if reference_path is not None:
                    reference = parse_alignment(reference_path)
                    prediction = match_prediction_names(prediction, reference)
                    quality = score(prediction, reference)
                if sci_config is not None:
                    sci_quality = calculate_sci(prediction, sci_config)
                    quality = {**(quality or {}), **sci_quality}
            except Exception as error:  # preserve failed scoring as evidence
                score_error = f"{type(error).__name__}: {error}"

        result = {
            "schema_version": 1,
            "method": "LinearTurboFold",
            "condition": condition_id,
            "ribosum": condition.get("ribosum", "not_applicable"),
            "dataset": dataset_id,
            "repetition": repetition,
            "status": "timeout" if timed_out else
                      ("success" if return_code == 0 else "failed"),
            "return_code": return_code,
            "timed_out": timed_out,
            "started_utc": started.isoformat(),
            "harness_wall_seconds": harness_seconds,
            "command": command,
            "command_sha256": digest,
            "input": str(input_path),
            "input_sha256": sha256_file(input_path),
            "reference": str(reference_path) if reference_path else None,
            "resource": load_resource(resource_path),
            "metrics_validation": {
                "valid": True,
                "not_applicable": True,
                "reason": "LinearTurboFold has no DAFS internal metrics trace",
            },
            "quality": quality,
            "quality_complete": (None if reference_path is None else
                                 bool(quality is not None and
                                      not quality["missing_predicted_sequences"] and
                                      not quality["extra_predicted_sequences"])),
            "score_error": score_error,
            "artifacts": {
                "alignment": str(alignment_path),
                "structures": str(run_dir),
                "stdout": str(stdout_path),
                "stderr": str(stderr_path),
                "resource": str(resource_path),
            },
        }
        atomic_json(result_path, result)
        print(f"{result['status'].upper()} {condition_id}/{dataset_id}/rep-{repetition:02d}")

    if args.jobs == 1:
        for job in jobs:
            run_job(job)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [executor.submit(run_job, job) for job in jobs]
            for future in concurrent.futures.as_completed(futures):
                future.result()

    if not args.dry_run:
        write_summaries(output_dir)
        print(f"Summary: {output_dir / 'summary.csv'}")
        failures = []
        for path in sorted((output_dir / "runs").glob("*/*/rep-*/result.json")):
            result = json.loads(path.read_text(encoding="utf-8"))
            if result.get("status") != "success":
                failures.append(f"{path}: status={result.get('status')}")
            elif result.get("score_error") and (
                    result.get("reference") or sci_config is not None):
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
