#!/usr/bin/env python3
"""Run external RNA aligners through one reproducible benchmark harness.

The runner deliberately keeps the command line for every method in JSON.  A
method may write an alignment to stdout (MAFFT, RAF, CentroidAlign) or to a
tool-specific result path (LocARNA).  CLUSTAL, Stockholm, and FASTA outputs
are normalized to ``prediction.fa`` before the common DAFS scorer is called.
For tools that do not emit a consensus structure, ``consensus_structure`` can
be set to ``rnaalifold`` to make the structure-prediction step explicit.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import platform
import random
import shlex
import shutil
import sys
import time
from collections import OrderedDict
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from run import (atomic_json, command_digest, execute, git_value,
                 load_resource, resolve_path, resolve_sci_config,
                 resolve_time_binary, sha256_file, validate_id,
                 write_summaries)
from score import (calculate_sci, match_prediction_names, parse_alignment,
                   score, with_consensus_structure)


def load_datasets(config: dict[str, Any], config_base: Path) -> list[dict[str, Any]]:
    datasets = [dict(dataset) for dataset in config.get("datasets", [])]
    for manifest_name in config.get("dataset_manifests", []):
        manifest_path = resolve_path(str(manifest_name), config_base)
        manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest_base = manifest_path.parent
        for dataset in manifest_data["datasets"]:
            resolved = dict(dataset)
            resolved["input"] = str(resolve_path(str(dataset["input"]), manifest_base))
            if dataset.get("reference"):
                resolved["reference"] = str(
                    resolve_path(str(dataset["reference"]), manifest_base))
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


def resolve_executable(value: str | Path, base: Path) -> Path:
    expanded = os.path.expandvars(str(value)).strip()
    candidate = Path(expanded).expanduser()
    if candidate.is_absolute() or "/" in expanded:
        resolved = candidate if candidate.is_absolute() else (base / candidate)
        executable = resolved.absolute()
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise FileNotFoundError(f"competitor executable not executable: {executable}")
        # Keep the invocation name: MAFFT dispatches Q-INS-i/X-INS-i/G-INS-i
        # through symlink aliases, so resolving the symlink changes the method.
        return executable
    found = shutil.which(expanded)
    if found is None:
        raise FileNotFoundError(
            f"competitor executable {expanded!r} was not found on PATH")
    return Path(found).absolute()


def expand_argument(value: Any, variables: dict[str, str]) -> str:
    if not isinstance(value, str):
        value = str(value)
    return os.path.expandvars(value).format_map(variables)


def method_binary(spec: dict[str, Any], config_base: Path) -> Path:
    value = spec.get("binary")
    if not value:
        raise ValueError(f"method {spec.get('id')!r} lacks binary")
    return resolve_executable(str(value), config_base)


def method_command(spec: dict[str, Any], binary: Path, input_path: Path,
                   run_dir: Path) -> list[str]:
    variables = {
        "binary": str(binary),
        "binary_dir": str(binary.parent),
        "input": str(input_path),
        "run_dir": str(run_dir),
        "input_stem": input_path.stem,
    }
    argv = spec.get("argv", ["{binary}", "{input}"])
    if not isinstance(argv, list) or not argv:
        raise ValueError(f"method {spec.get('id')!r} argv must be a non-empty list")
    return [expand_argument(value, variables) for value in argv]


def unsupported_input_symbols(input_path: Path,
                              supported_alphabet: str | list[str] | None
                              ) -> list[str]:
    """Return sequence symbols outside a method's explicitly supported alphabet."""
    if supported_alphabet is None:
        return []
    allowed = ({symbol.upper() for symbol in supported_alphabet}
               if isinstance(supported_alphabet, str) else
               {str(symbol).upper() for symbol in supported_alphabet})
    observed: set[str] = set()
    for line in input_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(">"):
            continue
        observed.update(character.upper() for character in stripped
                        if not character.isspace())
    return sorted(observed - allowed)


def output_source(spec: dict[str, Any], run_dir: Path,
                  stdout_path: Path) -> Path:
    mode = str(spec.get("output_mode", "stdout"))
    if mode == "stdout":
        return stdout_path
    if mode != "path":
        raise ValueError(
            f"method {spec.get('id')!r} output_mode must be stdout or path")
    variables = {"run_dir": str(run_dir)}
    output_path = spec.get("output_path")
    candidates: list[Path] = []
    if output_path:
        candidates.append(Path(expand_argument(output_path, variables)))
    for pattern in spec.get("output_glob", []):
        candidates.extend(Path(run_dir).glob(expand_argument(pattern, variables)))
    existing = [path if path.is_absolute() else path.resolve()
                for path in candidates if path.is_file()]
    if not existing:
        detail = output_path or spec.get("output_glob")
        raise FileNotFoundError(
            f"method {spec.get('id')!r} did not produce configured output {detail!r}")
    return existing[0]


def adapt_prediction(spec: dict[str, Any], prediction: Any) -> Any:
    """Apply explicit tool-output adapters before scoring."""
    dropped = set(map(str, spec.get("drop_sequence_names", [])))
    if not dropped:
        return prediction
    missing = dropped - set(prediction.sequences)
    if missing:
        raise ValueError(
            f"configured output-only records were not present: {sorted(missing)}")
    sequences = OrderedDict(
        (name, sequence) for name, sequence in prediction.sequences.items()
        if name not in dropped)
    structures = (None if prediction.structures is None else OrderedDict(
        (name, structure) for name, structure in prediction.structures.items()
        if name not in dropped))
    if not sequences:
        raise ValueError("output adapter removed every sequence")
    return replace(prediction, sequences=sequences, structures=structures)


def resolve_environment(config: dict[str, Any], config_base: Path,
                        sci_config: Any = None) -> dict[str, str]:
    """Resolve configured path-like environment values against the config."""
    environment = os.environ.copy()
    path_keys = {"MAFFT_BINARIES"}
    configured = config.get("environment", {})
    for raw_key, raw_value in configured.items():
        key = str(raw_key)
        value = os.path.expandvars(str(raw_value))
        if key.endswith(("_DIR", "_PATH")) or key in path_keys:
            candidate = Path(value).expanduser()
            if not candidate.is_absolute():
                value = str(resolve_path(value, config_base))
        environment[key] = value
    environment.setdefault("OMP_NUM_THREADS", "1")
    if sci_config is not None:
        vienna_bin = str(Path(sci_config.rnafold).parent)
        current_path = environment.get("PATH", "")
        environment["PATH"] = (vienna_bin + os.pathsep + current_path
                                if current_path else vienna_bin)
    return environment


def write_fasta_alignment(alignment: Any, destination: Path) -> None:
    with destination.open("w", encoding="utf-8") as stream:
        if alignment.structure is not None:
            stream.write(f">SS_cons\n{alignment.structure}\n")
        for name, sequence in alignment.sequences.items():
            stream.write(f">{name}\n{sequence}\n")


def clear_case_outputs(run_dir: Path) -> None:
    if not run_dir.exists():
        return
    for path in (run_dir / "prediction.fa", run_dir / "prediction.raw",
                 run_dir / "result.json", run_dir / "stdout.log",
                 run_dir / "stderr.log", run_dir / "resource.json"):
        if path.exists():
            path.unlink()
    tool_dir = run_dir / "tool"
    if tool_dir.exists():
        shutil.rmtree(tool_dir)


def run_one(spec: dict[str, Any], dataset: dict[str, Any], repetition: int,
            config: dict[str, Any], config_path: Path, config_base: Path,
            output_dir: Path, sci_config: Any, environment: dict[str, str],
            time_binary: Path | None,
            args: argparse.Namespace) -> None:
    method_id = validate_id(str(spec["id"]), "method")
    dataset_id = validate_id(str(dataset["id"]), "dataset")
    input_path = resolve_path(str(dataset["input"]), config_base)
    reference_path = (resolve_path(str(dataset["reference"]), config_base)
                      if dataset.get("reference") else None)
    if not input_path.is_file():
        raise FileNotFoundError(f"input not found: {input_path}")
    if reference_path is not None and not reference_path.is_file():
        raise FileNotFoundError(f"reference not found: {reference_path}")

    run_dir = output_dir / "runs" / method_id / dataset_id / f"rep-{repetition:02d}"
    prediction_path = run_dir / "prediction.fa"
    raw_path = run_dir / "prediction.raw"
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    resource_path = run_dir / "resource.json"
    result_path = run_dir / "result.json"
    binary = method_binary(spec, config_base)
    command = method_command(spec, binary, input_path, run_dir)
    digest = command_digest(command, sha256_file(input_path), {
        "runner": "run_competitor.py",
        "method": spec,
        "sci": (None if sci_config is None else {
            "rnaalifold": str(sci_config.rnaalifold),
            "rnafold": str(sci_config.rnafold),
            "version": sci_config.expected_version,
        }),
    })
    incompatible_symbols = unsupported_input_symbols(
        input_path, spec.get("input_alphabet"))

    if result_path.exists() and not args.rerun:
        previous = json.loads(result_path.read_text(encoding="utf-8"))
        completed_status = previous.get("status") in {
            "success", "skipped_incompatible_input"}
        has_artifact = (previous.get("status") == "skipped_incompatible_input"
                        or prediction_path.is_file())
        failed = not completed_status
        if (previous.get("command_sha256") == digest and has_artifact
                and not (args.rerun_failed and failed)):
            print(f"SKIP {method_id}/{dataset_id}/rep-{repetition:02d}")
            return

    if incompatible_symbols:
        reason = ("unsupported input symbols: "
                  + ", ".join(incompatible_symbols))
        print(f"SKIP-INCOMPATIBLE {method_id}/{dataset_id}: {reason}", flush=True)
        if args.dry_run:
            return
        run_dir.mkdir(parents=True, exist_ok=True)
        clear_case_outputs(run_dir)
        result = {
            "schema_version": 1,
            "method": str(spec.get("name", spec["id"])),
            "condition": method_id,
            "ribosum": "not_applicable",
            "dataset": dataset_id,
            "repetition": repetition,
            "status": "skipped_incompatible_input",
            "return_code": None,
            "timed_out": False,
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "harness_wall_seconds": 0.0,
            "command": command,
            "command_sha256": digest,
            "binary": str(binary),
            "binary_sha256": sha256_file(binary),
            "config_path": str(config_path),
            "config_sha256": sha256_file(config_path),
            "input": str(input_path),
            "input_sha256": sha256_file(input_path),
            "reference": str(reference_path) if reference_path else None,
            "reference_sha256": (sha256_file(reference_path)
                                 if reference_path is not None else None),
            "resource": {},
            "metrics_validation": {
                "valid": False,
                "kind": "external_competitor",
                "reason": "unsupported_input_alphabet",
                "unsupported_symbols": incompatible_symbols,
            },
            "quality": None,
            "quality_complete": None,
            "score_error": None,
            "skip_reason": reason,
            "artifacts": {
                "prediction": None,
                "raw_alignment": None,
                "stdout": None,
                "stderr": None,
                "resource": None,
            },
        }
        atomic_json(result_path, result)
        return

    print(f"RUN  {shlex.join(command)}", flush=True)
    if args.dry_run:
        return
    run_dir.mkdir(parents=True, exist_ok=True)
    clear_case_outputs(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    # All configured tools write under this directory; creating it here makes
    # LocARNA's --tgtdir and similar options independent of current cwd.
    (run_dir / "tool").mkdir(parents=True, exist_ok=True)
    time_format = ("{\"wall_seconds\":%e,\"user_seconds\":%U,"
                   "\"system_seconds\":%S,\"max_rss_kb\":%M,"
                   "\"exit_status\":%x}")
    wrapped = [str(time_binary), "-f", time_format, "-o", str(resource_path),
               "--", *command]
    started = datetime.now(timezone.utc)
    monotonic_start = time.monotonic()
    return_code, timed_out = execute(
        wrapped, stdout_path, stderr_path,
        float(config.get("timeout_seconds", 3600)),
        environment,
        cwd=(resolve_path(str(config["cwd"]), config_base)
             if config.get("cwd") else None))
    harness_seconds = time.monotonic() - monotonic_start
    quality = None
    score_error = None
    raw_source = None
    try:
        if return_code != 0:
            raise RuntimeError(f"competitor exited with status {return_code}")
        raw_source = output_source(spec, run_dir, stdout_path)
        shutil.copyfile(raw_source, raw_path)
        prediction = adapt_prediction(spec, parse_alignment(raw_source))
        if spec.get("consensus_structure") == "rnaalifold":
            if sci_config is None:
                raise ValueError(
                    "consensus_structure=rnaalifold requires a version-pinned sci configuration")
            prediction = with_consensus_structure(prediction, sci_config)
        write_fasta_alignment(prediction, prediction_path)
        reference = None
        if reference_path is not None:
            reference = parse_alignment(reference_path)
            prediction = match_prediction_names(prediction, reference)
            quality = score(prediction, reference)
        if sci_config is not None:
            quality = {**(quality or {}), **calculate_sci(prediction, sci_config)}
    except Exception as error:  # retain a complete failure record
        score_error = f"{type(error).__name__}: {error}"

    result = {
        "schema_version": 1,
        "method": str(spec.get("name", spec["id"])),
        "condition": method_id,
        "ribosum": "not_applicable",
        "dataset": dataset_id,
        "repetition": repetition,
        "status": ("timeout" if timed_out else
                   ("success" if return_code == 0 and score_error is None else "failed")),
        "return_code": return_code,
        "timed_out": timed_out,
        "started_utc": started.isoformat(),
        "harness_wall_seconds": harness_seconds,
        "command": command,
        "command_sha256": digest,
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "config_path": str(config_path),
        "config_sha256": sha256_file(config_path),
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "reference": str(reference_path) if reference_path else None,
        "reference_sha256": (sha256_file(reference_path)
                             if reference_path is not None else None),
        "resource": load_resource(resource_path),
        "metrics_validation": {
            "valid": return_code == 0 and score_error is None,
            "kind": "external_competitor",
        },
        "quality": quality,
        "quality_complete": (None if reference_path is None else
                              bool(quality is not None and
                                   not quality.get("missing_predicted_sequences") and
                                   not quality.get("extra_predicted_sequences"))),
        "score_error": score_error,
        "artifacts": {
            "prediction": str(prediction_path),
            "raw_alignment": str(raw_path) if raw_source is not None else None,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            "resource": str(resource_path),
        },
    }
    atomic_json(result_path, result)
    print(f"{result['status'].upper()} {method_id}/{dataset_id}/rep-{repetition:02d}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--method", action="append",
                        help="run only the named method id; may be repeated")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--rerun-failed", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--sci-bin-dir", type=Path)
    parser.add_argument("--sci-version")
    args = parser.parse_args()
    if args.jobs <= 0:
        parser.error("--jobs must be positive")

    config_path = args.config.resolve()
    config_base = config_path.parent
    config = json.loads(config_path.read_text(encoding="utf-8"))
    methods = [dict(spec) for spec in config.get("methods", [])]
    if not methods:
        raise ValueError("config.methods must contain at least one method")
    selected = set(args.method or [str(spec["id"]) for spec in methods])
    methods = [spec for spec in methods if str(spec["id"]) in selected]
    if not methods:
        raise ValueError("--method selected no configured method")
    for spec in methods:
        validate_id(str(spec["id"]), "method")

    sci_value = config.get("sci")
    if args.sci_bin_dir or args.sci_version:
        sci_value = dict(sci_value or {})
        if args.sci_bin_dir:
            sci_value["bin_dir"] = str(args.sci_bin_dir)
        if args.sci_version:
            sci_value["version"] = args.sci_version
    sci_config = resolve_sci_config(sci_value, config_base)
    repository = Path(__file__).resolve().parents[1]
    output_dir = (args.output_dir.resolve() if args.output_dir else
                  resolve_path(str(config.get("output_dir", "results-competitors")),
                               config_base))
    binaries = {str(spec["id"]): method_binary(spec, config_base)
                for spec in methods}
    time_binary = (None if args.dry_run else resolve_time_binary(
        config.get("time_binary", os.environ.get("TIME_BINARY", "time"))))
    datasets = load_datasets(config, config_base)
    repetitions = int(config.get("repetitions", 1))
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    environment = resolve_environment(config, config_base, sci_config)
    jobs = [(spec, dataset, repetition)
            for spec in methods for dataset in datasets
            for repetition in range(1, repetitions + 1)]
    random.Random(int(config.get("random_seed", 42))).shuffle(jobs)

    if not args.dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "runs").mkdir(exist_ok=True)
        manifest = {
            "schema_version": 1,
            "method_group": config.get("method_group", "external_competitors"),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "config": config,
            "config_path": str(config_path),
            "config_sha256": sha256_file(config_path),
            "binaries": {method_id: {
                "path": str(binary), "sha256": sha256_file(binary)}
                for method_id, binary in binaries.items()},
            "git_commit": git_value(repository, "rev-parse", "HEAD"),
            "git_status": git_value(repository, "status", "--short"),
            "platform": platform.platform(),
            "uname": list(platform.uname()),
            "python": sys.version,
            "environment": {key: environment[key]
                            for key in sorted(config.get("environment", {}))},
            "sci": (None if sci_config is None else {
                "rnaalifold": str(sci_config.rnaalifold),
                "rnafold": str(sci_config.rnafold),
                "expected_version": sci_config.expected_version,
                "rnaalifold_sha256": sha256_file(Path(sci_config.rnaalifold)),
                "rnafold_sha256": sha256_file(Path(sci_config.rnafold)),
            }),
        }
        atomic_json(output_dir / "manifest.json", manifest)

    def run_job(job: tuple[dict[str, Any], dict[str, Any], int]) -> None:
        run_one(job[0], job[1], job[2], config, config_path, config_base,
                output_dir, sci_config, environment, time_binary, args)

    if args.jobs == 1:
        for job in jobs:
            run_job(job)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [executor.submit(run_job, job) for job in jobs]
            for future in concurrent.futures.as_completed(futures):
                future.result()

    if args.dry_run:
        return 0
    write_summaries(output_dir)
    print(f"Summary: {output_dir / 'summary.csv'}")
    failures = []
    for path in sorted((output_dir / "runs").glob("*/*/rep-*/result.json")):
        result = json.loads(path.read_text(encoding="utf-8"))
        if result.get("status") not in {
                "success", "skipped_incompatible_input"}:
            failures.append(f"{path}: status={result.get('status')}")
    if failures and not config.get("allow_failures", False):
        for failure in failures:
            print(f"ERROR {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
