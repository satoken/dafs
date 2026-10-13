#!/usr/bin/env python3
"""Prepare and compare the round4 LinearPartition/LinearAlifold evidence.

The script only writes evaluation configs/reports and temporary probe binaries.
It never builds the production tree.  ``prepare`` is safe to run before a
candidate exists.  ``probe`` and ``compare-runs`` should be run after the
candidate build has been frozen and timing work has been coordinated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shlex
import struct
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from statistics import fmean
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"
CONFIG_DIR = BENCHMARKS / "configs"
PILOT_SOURCE = CONFIG_DIR / "murlet-dataset1-linearalifold-pilot-iter100-beam100.json"
FULL_SOURCE = CONFIG_DIR / "murlet-dataset1-linfold-optimization-iter100-beam100.json"
BASELINE_BINARY = Path("/tmp/dafs-luna-round4-baseline-freeze/dafs")
CANDIDATE_BINARY = Path("/tmp/dafs-luna-round4-candidate/src/dafs")
BASELINE_OBJECT_DIR = Path("/tmp/dafs-luna-round4-baseline-freeze/objects")
PROBE_SOURCE = BENCHMARKS / "linfold_round4_profile_probe.cpp"
VIENNA_BIN_DIR = Path("/home/sato-lab.org/satoken/.local/app/viennarna-2.7.2-arm64/bin")

PILOT_IDS = [
    "murlet-dataset1-RF00031-0",
    "murlet-dataset1-RF00003-0",
    "murlet-dataset1-RF00017-0",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False,
                   allow_nan=False) + "\n",
        encoding="utf-8",
    )


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def absolute(path: Path) -> Path:
    return path.expanduser().resolve()


def pair_config(source_path: Path, binary: Path, output_dir: Path,
                scope: str, dataset_ids: list[str] | None) -> dict[str, Any]:
    config = load_json(source_path)
    config["binary"] = str(absolute(binary))
    config["output_dir"] = str(absolute(output_dir))
    config["evaluation_scope"] = scope
    config["evaluation_source_config"] = str(absolute(source_path))
    config["repetitions"] = 1
    config.setdefault("environment", {})["OMP_NUM_THREADS"] = "1"
    if dataset_ids is not None:
        config["dataset_filter"] = {
            "collections": ["dataset1"],
            "ids": list(dataset_ids),
        }
    return config


def prepare(args: argparse.Namespace) -> int:
    baseline = absolute(args.baseline_binary)
    candidate = absolute(args.candidate_binary)
    pilot_source = absolute(args.pilot_source)
    full_source = absolute(args.full_source)
    specs = [
        ("pilot", pilot_source, PILOT_IDS),
        ("full", full_source, None),
    ]
    for scope, source, dataset_ids in specs:
        for label, binary in (("baseline", baseline), ("candidate", candidate)):
            config = pair_config(
                source,
                binary,
                Path(args.results_root) / f"{scope}-{label}",
                f"round4-{scope}",
                dataset_ids,
            )
            path = CONFIG_DIR / (
                f"murlet-dataset1-linfold-optimization-round4-"
                f"{scope}-iter100-beam100-{label}.json"
            )
            write_json(path, config)
            print(f"{label}_{scope}_config={path}")
            print(f"{label}_{scope}_config_sha256={sha256(path)}")
    print("pilot_expected_runs=24_per_variant")
    print("full_expected_runs=340_per_variant")
    print("run commands (not executed):")
    for scope in ("pilot", "full"):
        for label in ("baseline", "candidate"):
            path = CONFIG_DIR / (
                f"murlet-dataset1-linfold-optimization-round4-"
                f"{scope}-iter100-beam100-{label}.json"
            )
            print(f"  python3 benchmarks/run.py --config {shlex.quote(str(path))} --jobs 1")
    print("interleaved timing note: keep these correctness runs separate from any timing agent")
    return 0


def result_key(result: dict[str, Any]) -> tuple[str, str, int] | None:
    try:
        return (str(result["condition"]), str(result["dataset"]),
                int(result["repetition"]))
    except (KeyError, TypeError, ValueError):
        return None


def discover_results(root: Path) -> tuple[dict[tuple[str, str, int], tuple[Path, dict[str, Any]]], list[str]]:
    found: dict[tuple[str, str, int], tuple[Path, dict[str, Any]]] = {}
    errors: list[str] = []
    for path in sorted(root.glob("runs/*/*/rep-*/result.json")):
        try:
            result = load_json(path)
        except (OSError, json.JSONDecodeError) as error:
            errors.append(f"{path}: {error}")
            continue
        key = result_key(result)
        if key is None:
            errors.append(f"{path}: missing result identity")
        elif key in found:
            errors.append(f"{path}: duplicate identity {key}")
        else:
            found[key] = (path, result)
    return found, errors


TIME_KEYS = {
    "seconds", "dd_seconds", "wall_seconds", "user_seconds",
    "system_seconds", "cpu_seconds", "max_rss_kb", "harness_wall_seconds",
    "started_utc",
}
PATH_KEYS = {"input", "binary", "binary_path", "executable", "executable_path"}


def normalized(value: Any, *, trace: bool = False) -> Any:
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key in TIME_KEYS:
                continue
            if trace and key in PATH_KEYS:
                continue
            result[key] = normalized(item, trace=trace)
        return result
    if isinstance(value, list):
        return [normalized(item, trace=trace) for item in value]
    return value


def metric_trace(path: Path) -> tuple[list[Any], list[str]]:
    events: list[Any] = []
    errors: list[str] = []
    if not path.is_file():
        return events, [f"missing metrics: {path}"]
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        try:
            event = json.loads(line)
        except json.JSONDecodeError as error:
            errors.append(f"{path}:{line_number}: {error}")
            continue
        events.append(normalized(event, trace=True))
    return events, errors


def artifact_path(result: dict[str, Any], name: str) -> Path | None:
    value = (result.get("artifacts") or {}).get(name)
    return Path(value) if value else None


def compare_runs(args: argparse.Namespace) -> int:
    baseline_root = absolute(args.baseline_results)
    candidate_root = absolute(args.candidate_results)
    baseline, baseline_errors = discover_results(baseline_root)
    candidate, candidate_errors = discover_results(candidate_root)
    errors = [*baseline_errors, *candidate_errors]
    keys = sorted(set(baseline) | set(candidate))
    missing = [key for key in keys if key not in baseline or key not in candidate]
    records: list[dict[str, Any]] = []
    for key in keys:
        if key not in baseline or key not in candidate:
            continue
        base_path, base = baseline[key]
        cand_path, cand = candidate[key]
        record: dict[str, Any] = {
            "key": f"{key[0]}/{key[1]}/rep-{key[2]:02d}",
            "baseline_status": base.get("status"),
            "candidate_status": cand.get("status"),
        }
        base_prediction = artifact_path(base, "prediction")
        cand_prediction = artifact_path(cand, "prediction")
        record["prediction_equal"] = bool(
            base_prediction and cand_prediction and
            base_prediction.is_file() and cand_prediction.is_file() and
            base_prediction.read_bytes() == cand_prediction.read_bytes()
        )
        record["quality_equal"] = base.get("quality") == cand.get("quality")
        base_summary = normalized((base.get("metrics_validation") or {}).get("run_summary"))
        cand_summary = normalized((cand.get("metrics_validation") or {}).get("run_summary"))
        record["summary_non_time_equal"] = base_summary == cand_summary
        base_metrics = artifact_path(base, "metrics")
        cand_metrics = artifact_path(cand, "metrics")
        base_trace, base_trace_errors = metric_trace(base_metrics) if base_metrics else ([], [f"{base_path}: no metrics path"])
        cand_trace, cand_trace_errors = metric_trace(cand_metrics) if cand_metrics else ([], [f"{cand_path}: no metrics path"])
        errors.extend(base_trace_errors)
        errors.extend(cand_trace_errors)
        record["trace_equal"] = base_trace == cand_trace
        record["metrics_valid"] = bool(
            (base.get("metrics_validation") or {}).get("valid") and
            (cand.get("metrics_validation") or {}).get("valid")
        )
        for label, result in (("baseline", base), ("candidate", cand)):
            resource = result.get("resource") or {}
            for field in ("wall_seconds", "user_seconds", "system_seconds",
                          "max_rss_kb"):
                try:
                    record[f"{label}_{field}"] = float(resource[field])
                except (KeyError, TypeError, ValueError):
                    record[f"{label}_{field}"] = None
        records.append(record)

    def count(name: str, value: Any = True) -> int:
        return sum(record.get(name) == value for record in records)

    timing: dict[str, Any] = {}
    for label in ("baseline", "candidate"):
        for field in ("wall_seconds", "user_seconds", "system_seconds", "max_rss_kb"):
            values = [record[f"{label}_{field}"] for record in records
                      if record[f"{label}_{field}"] is not None]
            timing[f"{label}_{field}_mean"] = fmean(values) if values else None
            timing[f"{label}_{field}_max"] = max(values) if values else None

    report = {
        "schema_version": 1,
        "driver": str(Path(__file__).resolve()),
        "baseline_results": str(baseline_root),
        "candidate_results": str(candidate_root),
        "baseline_manifest": load_json(baseline_root / "manifest.json") if (baseline_root / "manifest.json").is_file() else None,
        "candidate_manifest": load_json(candidate_root / "manifest.json") if (candidate_root / "manifest.json").is_file() else None,
        "expected_runs": len(keys),
        "paired_runs": len(records),
        "missing_runs": [f"{key[0]}/{key[1]}/rep-{key[2]:02d}" for key in missing],
        "errors": errors,
        "counts": {
            "prediction_equal": count("prediction_equal"),
            "quality_equal": count("quality_equal"),
            "summary_non_time_equal": count("summary_non_time_equal"),
            "trace_equal": count("trace_equal"),
            "metrics_valid_both": count("metrics_valid"),
            "baseline_success": sum(record["baseline_status"] == "success" for record in records),
            "candidate_success": sum(record["candidate_status"] == "success" for record in records),
        },
        "timing_and_rss": timing,
        "records": records,
    }
    report_path = absolute(args.report)
    write_json(report_path, report)
    print(f"paired={len(records)} missing={len(missing)} errors={len(errors)}")
    print(json.dumps(report["counts"], sort_keys=True))
    print(f"report={report_path}")
    fatal = bool(errors or missing or any(
        not record[key]
        for record in records
        for key in ("prediction_equal", "quality_equal", "summary_non_time_equal", "trace_equal", "metrics_valid")
    ))
    return 1 if fatal else 0


def object_for(build: Path, name: str) -> Path:
    candidates = [
        build / "src" / "CMakeFiles" / "linear_profile_oracle_test.dir" / name,
        build / "src" / "CMakeFiles" / "dafs.dir" / name,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"could not find {name} under {build}")


def cmake_source_root(build: Path) -> Path:
    """Return the source tree used for a candidate build.

    Attribution builds may come from temporary source copies (for example a
    BPP-only copy with the frozen round3 profile header).  Compiling the probe
    against the live checkout in that case can silently give the wrapper a
    different Options layout.  Keep the probe headers paired with the object
    files that it links.
    """
    cache = build / "CMakeCache.txt"
    if cache.is_file():
        prefix = "CMAKE_HOME_DIRECTORY:INTERNAL="
        for line in cache.read_text(encoding="utf-8").splitlines():
            if line.startswith(prefix):
                return Path(line[len(prefix):]).resolve()
    return ROOT


def cmake_linfold_defines(build: Path) -> list[str]:
    """Extract the dense LinFold feature flags from an isolated build."""
    flags = build / "src" / "CMakeFiles" / "dafs.dir" / "flags.make"
    if not flags.is_file():
        return []
    text = flags.read_text(encoding="utf-8")
    defines: list[str] = []
    for name in ("BPP", "N_ALPHA", "N_BETA"):
        marker = f"DAFS_LINFOLD_DENSE_{name}="
        for token in text.split():
            if token.startswith("-D" + marker):
                defines.append(token)
                break
    return defines


def probe_command(variant: str, build: Path, output: Path) -> tuple[list[str], Path]:
    cxx = Path("/usr/bin/c++")
    if not cxx.is_file():
        cxx = Path("c++")
    work = Path("/tmp/dafs-luna-round4-profile-probe") / variant
    work.mkdir(parents=True, exist_ok=True)
    probe_object = work / "probe.o"
    executable = work / "profile_probe"
    source_root = (BASELINE_SOURCE if variant == "baseline"
                   else cmake_source_root(build))
    vienna_include = Path("/home/sato-lab.org/satoken/.local/app/viennarna-2.7.2-arm64/include")
    vienna_lib = Path("/home/sato-lab.org/satoken/.local/app/viennarna-2.7.2-arm64/lib")
    include_flags = [
        "-I", str(source_root / "src"),
        # The frozen source snapshot intentionally contains only the
        # LinFold/wrapper closure.  Fall back to the checkout for unrelated
        # public headers (sparse_matrix.h, fa.h, and their dependencies),
        # while keeping matching LinFold headers first.
        "-I", str(ROOT / "src"),
        "-I", str(ROOT / "spdlog/include"),
        "-I", str(ROOT / "cxxopts/include"),
        "-I", str(source_root / "src/linearfold/fold"),
        "-I", str(ROOT / "src/linearfold/fold"),
        "-I", str(vienna_include),
        "-I", str(vienna_include / "ViennaRNA"),
    ]
    if variant == "baseline":
        wrapper_object = BASELINE_OBJECT_DIR / "linfold_wrapper.cpp.o"
        ribosum_object = BASELINE_OBJECT_DIR / "ribosum.cpp.o"
        library = BASELINE_OBJECT_DIR / "liblinearfold.a"
    else:
        wrapper_object = object_for(build, "linfold_wrapper.cpp.o")
        ribosum_object = object_for(build, "ribosum.cpp.o")
        library = build / "src" / "linearfold" / "liblinearfold.a"
    if not wrapper_object.is_file() or not ribosum_object.is_file() or not library.is_file():
        raise FileNotFoundError(f"probe link inputs are incomplete for {variant}")
    compile_cmd = [str(cxx), "-O3", "-DNDEBUG", "-std=gnu++17",
                   *cmake_linfold_defines(build),
                   *include_flags, "-c", str(PROBE_SOURCE), "-o", str(probe_object)]
    link_cmd = [str(cxx), "-O3", "-DNDEBUG", "-fno-lto", "-Wl,-fno-lto",
                str(probe_object), str(wrapper_object), str(ribosum_object),
                "-o", str(executable), str(library), "-L", str(vienna_lib),
                "-lRNA", "-lgomp", "-lpthread", "-ldl"]
    subprocess.run(compile_cmd, check=True)
    subprocess.run(link_cmd, check=True)
    return [str(executable)], output


def run_probe(args: argparse.Namespace) -> int:
    candidate_build = absolute(args.candidate_build)
    for variant, build, output in (
        ("baseline", Path("/tmp/dafs-luna-round4-baseline-freeze"), absolute(args.baseline_output)),
        ("candidate", candidate_build, absolute(args.candidate_output)),
    ):
        command, output_path = probe_command(variant, build, output)
        print(f"probe_{variant}_command={shlex.join(command)}")
        with output_path.open("w", encoding="utf-8") as stream:
            subprocess.run(command, check=True, stdout=stream)
        print(f"probe_{variant}_output={output_path}")
    return 0


def float_from_bits(bits_text: str) -> float:
    return struct.unpack(">f", int(bits_text, 16).to_bytes(4, "big"))[0]


def load_probe(path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    records: dict[tuple[str, str], dict[str, Any]] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            value = json.loads(line)
            key = (str(value["run"]), str(value["case"]))
        except (json.JSONDecodeError, KeyError, TypeError) as error:
            raise ValueError(f"{path}:{line_number}: invalid probe row: {error}") from error
        if key in records:
            raise ValueError(f"{path}:{line_number}: duplicate probe row {key}")
        records[key] = value
    return records


def compare_probe_row(base: dict[str, Any], cand: dict[str, Any]) -> dict[str, Any]:
    base_pairs = {(int(pair["i"]), int(pair["j"])): pair for pair in base["pairs"]}
    cand_pairs = {(int(pair["i"]), int(pair["j"])): pair for pair in cand["pairs"]}
    common = set(base_pairs) & set(cand_pairs)
    max_abs = 0.0
    max_rel = 0.0
    changed_bits = 0
    for key in common:
        old = float_from_bits(base_pairs[key]["p_bits"])
        new = float_from_bits(cand_pairs[key]["p_bits"])
        delta = abs(old - new)
        max_abs = max(max_abs, delta)
        max_rel = max(max_rel, delta / max(abs(old), abs(new), 1e-30))
        changed_bits += base_pairs[key]["p_bits"] != cand_pairs[key]["p_bits"]
    return {
        "logz_bitwise_equal": base["logz_bits"] == cand["logz_bits"],
        "pair_keys_equal": set(base_pairs) == set(cand_pairs),
        "pair_bitwise_equal": changed_bits == 0 and set(base_pairs) == set(cand_pairs),
        "baseline_pair_count": len(base_pairs),
        "candidate_pair_count": len(cand_pairs),
        "common_pair_count": len(common),
        "pair_bit_changes": changed_bits,
        "max_abs_bpp_diff": max_abs,
        "max_rel_bpp_diff": max_rel,
    }


def compare_probes(args: argparse.Namespace) -> int:
    baseline = load_probe(absolute(args.baseline_output))
    candidate = load_probe(absolute(args.candidate_output))
    keys = sorted(set(baseline) | set(candidate))
    rows: list[dict[str, Any]] = []
    for key in keys:
        if key not in baseline or key not in candidate:
            rows.append({"run": key[0], "case": key[1], "missing": True})
            continue
        rows.append({"run": key[0], "case": key[1],
                     **compare_probe_row(baseline[key], candidate[key])})

    # Compare cutoff crossings using the complete cutoff=0 output, which keeps
    # all entries regardless of the reporting threshold.
    crossings: list[dict[str, Any]] = []
    for run_name, case_name in sorted({key for key in baseline if key[0] == "lpv_cutoff0"}):
        base_row = baseline[(run_name, case_name)]
        cand_row = candidate.get((run_name, case_name))
        if cand_row is None:
            continue
        base_map = {(int(pair["i"]), int(pair["j"])): float_from_bits(pair["p_bits"])
                    for pair in base_row["pairs"]}
        cand_map = {(int(pair["i"]), int(pair["j"])): float_from_bits(pair["p_bits"])
                    for pair in cand_row["pairs"]}
        for cutoff in (0.001, 0.01, 0.05, 0.1):
            keys_union = set(base_map) | set(cand_map)
            crossed = sorted(key for key in keys_union
                             if (base_map.get(key, 0.0) > cutoff) !=
                             (cand_map.get(key, 0.0) > cutoff))
            crossings.append({
                "case": case_name,
                "cutoff": cutoff,
                "baseline_above": sum(value > cutoff for value in base_map.values()),
                "candidate_above": sum(value > cutoff for value in cand_map.values()),
                "crossing_count": len(crossed),
                "first_crossings": [list(key) for key in crossed[:10]],
            })

    report = {
        "schema_version": 1,
        "driver": str(Path(__file__).resolve()),
        "baseline_output": str(absolute(args.baseline_output)),
        "candidate_output": str(absolute(args.candidate_output)),
        "baseline_binary_sha256": sha256(BASELINE_BINARY),
        "candidate_build": str(absolute(args.candidate_build)),
        "paired_rows": len([row for row in rows if not row.get("missing")]),
        "missing_rows": [row for row in rows if row.get("missing")],
        "rows": rows,
        "cutoff_crossings": crossings,
        "all_logz_bitwise_equal": all(row.get("logz_bitwise_equal") for row in rows if not row.get("missing")),
        "all_pair_keys_equal": all(row.get("pair_keys_equal") for row in rows if not row.get("missing")),
        "all_pair_bitwise_equal": all(row.get("pair_bitwise_equal") for row in rows if not row.get("missing")),
        "max_abs_bpp_diff": max((row.get("max_abs_bpp_diff", 0.0) for row in rows), default=0.0),
        "max_rel_bpp_diff": max((row.get("max_rel_bpp_diff", 0.0) for row in rows), default=0.0),
        "cutoff_crossing_total": sum(item["crossing_count"] for item in crossings),
    }
    report_path = absolute(args.report)
    write_json(report_path, report)
    print(json.dumps({
        "paired_rows": report["paired_rows"],
        "missing_rows": len(report["missing_rows"]),
        "all_logz_bitwise_equal": report["all_logz_bitwise_equal"],
        "all_pair_keys_equal": report["all_pair_keys_equal"],
        "all_pair_bitwise_equal": report["all_pair_bitwise_equal"],
        "max_abs_bpp_diff": report["max_abs_bpp_diff"],
        "max_rel_bpp_diff": report["max_rel_bpp_diff"],
        "cutoff_crossing_total": report["cutoff_crossing_total"],
    }, sort_keys=True))
    print(f"report={report_path}")
    return 1 if (
        report["missing_rows"] or
        not report["all_logz_bitwise_equal"] or
        not report["all_pair_keys_equal"] or
        not report["all_pair_bitwise_equal"]
    ) else 0


BASELINE_SOURCE = Path("/tmp/dafs-luna-round4-baseline-source")


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)

    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--baseline-binary", type=Path, default=BASELINE_BINARY)
    prepare_parser.add_argument("--candidate-binary", type=Path, default=CANDIDATE_BINARY)
    prepare_parser.add_argument("--pilot-source", type=Path, default=PILOT_SOURCE)
    prepare_parser.add_argument("--full-source", type=Path, default=FULL_SOURCE)
    prepare_parser.add_argument("--results-root", type=Path, default=Path("/tmp/dafs-luna-round4"))
    prepare_parser.set_defaults(handler=prepare)

    compare_parser = commands.add_parser("compare-runs")
    compare_parser.add_argument("--baseline-results", type=Path, required=True)
    compare_parser.add_argument("--candidate-results", type=Path, required=True)
    compare_parser.add_argument("--report", type=Path,
                                default=Path("/tmp/dafs-luna-round4-pilot-compare.json"))
    compare_parser.set_defaults(handler=compare_runs)

    probe_parser = commands.add_parser("probe")
    probe_parser.add_argument("--candidate-build", type=Path, required=True)
    probe_parser.add_argument("--baseline-output", type=Path,
                              default=Path("/tmp/dafs-luna-round4-profile-baseline.jsonl"))
    probe_parser.add_argument("--candidate-output", type=Path,
                              default=Path("/tmp/dafs-luna-round4-profile-candidate.jsonl"))
    probe_parser.set_defaults(handler=run_probe)

    probe_compare = commands.add_parser("compare-probes")
    probe_compare.add_argument("--candidate-build", type=Path, required=True)
    probe_compare.add_argument("--baseline-output", type=Path,
                               default=Path("/tmp/dafs-luna-round4-profile-baseline.jsonl"))
    probe_compare.add_argument("--candidate-output", type=Path,
                               default=Path("/tmp/dafs-luna-round4-profile-candidate.jsonl"))
    probe_compare.add_argument("--report", type=Path,
                               default=Path("/tmp/dafs-luna-round4-profile-compare.json"))
    probe_compare.set_defaults(handler=compare_probes)
    return root


def main() -> int:
    args = parser().parse_args()
    try:
        return int(args.handler(args))
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
