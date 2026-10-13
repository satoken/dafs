#!/usr/bin/env python3
"""Small round2-only harness for paired timing and sparse decoder scaling.

The full Murlet run and trace comparison stay in benchmark_luna_optimization.py
and run.py.  This file only owns the metrics-off serial pair and the isolated
Beam/LinearNussinov microbenchmark orchestration.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"
CONFIG_DIR = BENCHMARKS / "configs"
DEFAULT_CONFIG = CONFIG_DIR / "luna-round2-optimized.json"
DEFAULT_RESULTS = BENCHMARKS / "results-murlet-dataset1-luna-round2-arm64-v272"
BASELINE_BINARY = Path("/tmp/dafs-luna-round2-7MkkIm/build-arm64/src/dafs")
BASELINE_SHA256 = "49f9abb54daef26701ef0b396192b113c47ee851ddf6c86d23f5d687eee43432"
BACKUP_SOURCE = Path("/tmp/dafs-luna-round2-7MkkIm/src")
SCALING_SOURCE = BENCHMARKS / "luna_round2_scaling.cpp"
SERIAL_DATASET = "murlet-dataset1-RF00169-2"
SERIAL_CONDITION = "cl_linear_ath02_iter1000"
DEFAULT_LENGTHS = (250, 500, 1000, 2000, 4000)


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
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def ensure_fresh(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise RuntimeError(f"refusing to overwrite non-empty round2 path: {path}")
    path.mkdir(parents=True, exist_ok=True)


def read_config(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"config is not an object: {path}")
    return value


def serial_paths(config: dict[str, Any], config_path: Path) -> tuple[Path, list[str]]:
    manifests = config.get("dataset_manifests") or []
    if not manifests:
        raise ValueError("round2 config has no dataset manifest")
    manifest_path = Path(manifests[0]).expanduser()
    if not manifest_path.is_absolute():
        manifest_path = (config_path.parent / manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    dataset = next(item for item in manifest["datasets"]
                   if item["id"] == SERIAL_DATASET)
    input_path = (manifest_path.parent / dataset["input"]).resolve()
    condition = next(item for item in config["conditions"]
                     if item["id"] == SERIAL_CONDITION)
    return input_path, [str(argument) for argument in condition["args"]]


def timed_command(binary: Path, arguments: list[str], input_path: Path,
                  run_dir: Path, timeout_seconds: float,
                  environment: dict[str, str]) -> dict[str, Any]:
    run_dir.mkdir(parents=True, exist_ok=True)
    prediction = run_dir / "prediction.fa"
    stderr = run_dir / "stderr.log"
    resource = run_dir / "resource.json"
    command = [str(binary), *arguments, str(input_path)]
    time_binary = Path("/usr/bin/time")
    if not time_binary.is_file():
        found = shutil.which("time")
        if not found:
            raise RuntimeError("GNU time not found")
        time_binary = Path(found)
    wrapped = [str(time_binary), "-f",
               '{"wall_seconds":%e,"user_seconds":%U,'
               '"system_seconds":%S,"max_rss_kb":%M,"exit_status":%x}',
               "-o", str(resource), "--", *command]
    started = utc_now()
    try:
        with prediction.open("wb") as output, stderr.open("wb") as error:
            completed = subprocess.run(
                wrapped, cwd=ROOT, env=environment, stdout=output, stderr=error,
                timeout=timeout_seconds, check=False)
        return_code = completed.returncode
        timed_out = False
    except subprocess.TimeoutExpired:
        return_code = -9
        timed_out = True
    resource_value: dict[str, Any] = {}
    if resource.is_file():
        resource_value = json.loads(resource.read_text(encoding="utf-8"))
    result = {
        "schema_version": 1,
        "started_utc": started,
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "command": command,
        "command_sha256": hashlib.sha256(
            json.dumps(command, separators=(",", ":")).encode()).hexdigest(),
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "status": "timeout" if timed_out else
                  ("success" if return_code == 0 else "failed"),
        "return_code": return_code,
        "metrics_enabled": False,
        "resource": resource_value,
        "prediction_sha256": sha256_file(prediction) if prediction.is_file() else None,
        "artifacts": {"prediction": str(prediction), "stderr": str(stderr),
                      "resource": str(resource)},
    }
    write_json(run_dir / "result.json", result)
    return result


def run_serial(args: argparse.Namespace) -> int:
    config_path = Path(args.config).resolve()
    config = read_config(config_path)
    root = Path(args.results_root).resolve()
    serial_root = root / "serial-paired"
    ensure_fresh(serial_root)
    input_path, condition_args = serial_paths(config, config_path)
    baseline = Path(args.baseline_binary).resolve()
    optimized = Path(args.optimized_binary).resolve()
    if sha256_file(baseline) != BASELINE_SHA256:
        raise RuntimeError(f"baseline SHA256 mismatch: {baseline}")
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    rows: list[dict[str, Any]] = []
    repetitions = int(args.repetitions)
    for repetition in range(1, repetitions + 1):
        order = ([("baseline", baseline), ("optimized", optimized)]
                 if repetition % 2 else
                 [("optimized", optimized), ("baseline", baseline)])
        results: dict[str, dict[str, Any]] = {}
        for order_index, (side, binary) in enumerate(order, 1):
            run_dir = serial_root / side / f"rep-{repetition:02d}"
            result = timed_command(binary, condition_args, input_path, run_dir,
                                   float(args.timeout_seconds), environment)
            result["side"] = side
            result["repetition"] = repetition
            result["order_index"] = order_index
            results[side] = result
        old = results["baseline"]
        new = results["optimized"]
        old_resource = old.get("resource", {})
        new_resource = new.get("resource", {})
        old_wall = float(old_resource.get("wall_seconds", "nan"))
        new_wall = float(new_resource.get("wall_seconds", "nan"))
        old_cpu = (float(old_resource.get("user_seconds", 0)) +
                   float(old_resource.get("system_seconds", 0)))
        new_cpu = (float(new_resource.get("user_seconds", 0)) +
                   float(new_resource.get("system_seconds", 0)))
        row = {
            "repetition": repetition,
            "order": "baseline-optimized" if repetition % 2 else "optimized-baseline",
            "baseline_status": old["status"],
            "optimized_status": new["status"],
            "baseline_wall_seconds": old_wall,
            "optimized_wall_seconds": new_wall,
            "baseline_user_seconds": old_resource.get("user_seconds"),
            "optimized_user_seconds": new_resource.get("user_seconds"),
            "baseline_system_seconds": old_resource.get("system_seconds"),
            "optimized_system_seconds": new_resource.get("system_seconds"),
            "baseline_cpu_seconds": old_cpu,
            "optimized_cpu_seconds": new_cpu,
            "baseline_max_rss_kb": old_resource.get("max_rss_kb"),
            "optimized_max_rss_kb": new_resource.get("max_rss_kb"),
            "wall_speedup": old_wall / new_wall if new_wall > 0 else None,
            "cpu_speedup": old_cpu / new_cpu if new_cpu > 0 else None,
            "prediction_bytes_equal": old.get("prediction_sha256") ==
                                      new.get("prediction_sha256"),
        }
        rows.append(row)
    fieldnames = list(rows[0])
    with (serial_root / "paired.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    write_json(serial_root / "paired.json", {
        "schema_version": 1,
        "created_utc": utc_now(),
        "config": str(config_path),
        "dataset": SERIAL_DATASET,
        "condition": SERIAL_CONDITION,
        "repetitions": repetitions,
        "metrics_enabled": False,
        "rows": rows,
    })
    print(serial_root / "paired.csv")
    failed = any(
        row["baseline_status"] != "success" or
        row["optimized_status"] != "success" or
        not row["prediction_bytes_equal"]
        for row in rows
    )
    return 1 if failed else 0


def compile_one(command: list[str]) -> None:
    print("COMPILE", shlex.join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def compile_scaling(root: Path) -> dict[str, Any]:
    ensure_fresh(root)
    executable_dir = root / "executables"
    executable_dir.mkdir(parents=True)
    compiler = os.environ.get("CXX", "c++")
    common = [compiler, "-O3", "-DNDEBUG", "-std=c++17"]
    executables: dict[str, dict[str, Any]] = {}

    def build(name: str, side: str, macro: str, source: Path,
              includes: list[str], extra: list[str] = []) -> None:
        output = executable_dir / name
        command = common + [f"-D{macro}"]
        for include in includes:
            command += [f"-I{include}"]
        command += extra + [str(SCALING_SOURCE), str(source), "-o", str(output)]
        compile_one(command)
        executables[name] = {
            "side": side, "path": str(output), "sha256": sha256_file(output),
            "command": command,
        }

    workspace_src = ROOT / "src"
    backup_src = BACKUP_SOURCE
    build("baseline-beam", "baseline", "LUNA_ROUND2_BEAM",
          backup_src / "linearalign/BeamAlign.cpp",
          [str(backup_src / "linearalign"), str(workspace_src / "linearalign"),
           str(workspace_src)])
    build("optimized-beam", "optimized", "LUNA_ROUND2_BEAM",
          workspace_src / "linearalign/BeamAlign.cpp",
          [str(workspace_src / "linearalign"), str(workspace_src)])
    build("baseline-nuss", "baseline", "LUNA_ROUND2_NUSS",
          backup_src / "nussinov.cpp", [str(workspace_src)],
          ["-include", str(workspace_src / "nussinov.h")])
    build("optimized-nuss", "optimized", "LUNA_ROUND2_NUSS",
          workspace_src / "nussinov.cpp", [str(workspace_src)])
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "compiler": compiler,
        "workspace": str(ROOT),
        "backup_source": str(BACKUP_SOURCE),
        "scaling_source": str(SCALING_SOURCE),
        "scaling_source_sha256": sha256_file(SCALING_SOURCE),
        "executables": executables,
        "host": {"platform": platform.platform(), "uname": list(os.uname())},
    }
    write_json(root / "compile.json", manifest)
    return manifest


def read_one_scaling_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise RuntimeError(f"expected one scaling row in {path}, got {len(rows)}")
    return rows[0]


def run_scaling(root: Path, lengths: tuple[int, ...], repetitions: int,
                warmup: int) -> int:
    manifest = json.loads((root / "compile.json").read_text(encoding="utf-8"))
    result_root = root / "results"
    ensure_fresh(result_root)
    variants = (("beam", "baseline-beam", "optimized-beam"),
                ("nuss-certified", "baseline-nuss", "optimized-nuss"),
                ("nuss-uncertified", "baseline-nuss", "optimized-nuss"))
    rows: list[dict[str, Any]] = []
    invocation = 0
    for variant, baseline_name, optimized_name in variants:
        for length in lengths:
            for side, executable_name in (("baseline", baseline_name),
                                           ("optimized", optimized_name)):
                invocation += 1
                executable = Path(manifest["executables"][executable_name]["path"])
                row_path = result_root / f"{invocation:02d}-{side}-{variant}-{length}.csv"
                command = [str(executable), "--variant", variant,
                           "--lengths", str(length), "--beam", "100",
                           "--repetitions", str(repetitions), "--warmup", str(warmup),
                           "--seed", "42", "--output", str(row_path)]
                print("RUN", shlex.join(command), flush=True)
                subprocess.run(command, cwd=ROOT, check=True)
                row = read_one_scaling_row(row_path)
                row.update({"side": side, "executable": str(executable),
                            "executable_sha256": sha256_file(executable),
                            "invocation": invocation})
                rows.append(row)
    fieldnames = list(rows[0])
    with (result_root / "scaling_rows.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    by_key = {(row["variant"], row["length"], row["side"]): row for row in rows}
    comparisons: list[dict[str, Any]] = []
    for variant, _, _ in variants:
        for length in lengths:
            old = by_key[(variant, str(length), "baseline")]
            new = by_key[(variant, str(length), "optimized")]
            digest_fields = ("params_digest", "input_digest", "score_digest",
                             "structure_digest", "bounds_digest")
            digest_equal = {field: old[field] == new[field] for field in digest_fields}
            def number(row: dict[str, str], field: str) -> float | None:
                value = row.get(field, "")
                return None if value in ("", "nan", "NaN") else float(value)
            old_wall, new_wall = number(old, "wall_seconds"), number(new, "wall_seconds")
            old_cpu, new_cpu = number(old, "cpu_seconds"), number(new, "cpu_seconds")
            comparisons.append({
                "variant": variant, "length": length,
                "support_mode": old["support_mode"],
                "beam": int(old["beam"]),
                "p_nnz": int(old["p_nnz"]), "q_nnz": int(old["q_nnz"]),
                "raw_support": int(old["raw_support"]),
                "unique_support": int(old["unique_support"]),
                "digest_equal": digest_equal,
                "all_digest_equal": all(digest_equal.values()),
                "output_stable": old["output_stable"] == "true" and
                                 new["output_stable"] == "true",
                "baseline_score": number(old, "score"),
                "optimized_score": number(new, "score"),
                "baseline_wall_seconds": old_wall,
                "optimized_wall_seconds": new_wall,
                "baseline_cpu_seconds": old_cpu,
                "optimized_cpu_seconds": new_cpu,
                "baseline_max_rss_kb": int(old["max_rss_kb"]),
                "optimized_max_rss_kb": int(new["max_rss_kb"]),
                "wall_speedup": old_wall / new_wall if old_wall and new_wall else None,
                "cpu_speedup": old_cpu / new_cpu if old_cpu and new_cpu else None,
            })
    compare_fields = list(comparisons[0])
    with (result_root / "scaling_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=compare_fields)
        writer.writeheader()
        writer.writerows(comparisons)
    write_json(result_root / "scaling_comparison.json", {
        "schema_version": 1, "created_utc": utc_now(),
        "lengths": list(lengths), "beam": 100, "repetitions": repetitions,
        "warmup": warmup,
        "rows": comparisons,
    })
    print(result_root / "scaling_comparison.csv")
    failed = any(not row["all_digest_equal"] or not row["output_stable"]
                 for row in comparisons)
    return 1 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    serial = subparsers.add_parser("serial", help="metrics-off paired subset")
    serial.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    serial.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS)
    serial.add_argument("--baseline-binary", type=Path, default=BASELINE_BINARY)
    serial.add_argument("--optimized-binary", type=Path,
                        default=ROOT / "build-arm64/src/dafs")
    serial.add_argument("--repetitions", type=int, choices=(3, 5), default=5)
    serial.add_argument("--timeout-seconds", type=float, default=3600)
    serial.set_defaults(handler=run_serial)

    scaling = subparsers.add_parser("scaling", help="compile/run/compare decoder scale")
    scaling.add_argument("action", choices=("compile", "run", "all"))
    scaling.add_argument("--root", type=Path,
                         default=DEFAULT_RESULTS / "scaling")
    scaling.add_argument("--lengths", type=int, nargs="+", default=list(DEFAULT_LENGTHS))
    scaling.add_argument("--repetitions", type=int, default=3)
    scaling.add_argument("--warmup", type=int, default=1)
    scaling.set_defaults(handler=None)

    args = parser.parse_args()
    try:
        if args.command == "scaling":
            root = args.root.resolve()
            if args.action in ("compile", "all"):
                compile_scaling(root)
            if args.action in ("run", "all"):
                return run_scaling(root, tuple(args.lengths), args.repetitions,
                                   args.warmup)
            return 0
        return int(args.handler(args))
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
