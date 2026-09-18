#!/usr/bin/env python3
"""Create a compact, machine-readable bundle from a DAFS benchmark run."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import tarfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

NUMERIC_FIELDS = (
    "wall_seconds", "user_seconds", "system_seconds", "max_rss_kb",
    "sps", "alignment_ppv", "sensitivity", "ppv", "mcc", "cbp_f1", "sci",
    "sci_consensus_mfe", "sci_mean_single_mfe", "internal_seconds",
    "dd_iterations", "cbp_peak",
)


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    rows = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if line.strip():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError as error:
                raise ValueError(f"{path}:{line_number}: {error}") from error
    return rows


def numeric(values: Iterable[Any]) -> list[float]:
    return [float(value) for value in values
            if isinstance(value, (int, float)) and not isinstance(value, bool)]


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "runs": len(rows),
        "successes": sum(row.get("status") == "success" for row in rows),
        "timeouts": sum(row.get("status") == "timeout" for row in rows),
        "failures": sum(row.get("status") == "failed" for row in rows),
        "metrics_valid": sum(row.get("metrics_valid") is True for row in rows),
        "quality_complete": sum(row.get("quality_complete") is True for row in rows),
    }
    for field in NUMERIC_FIELDS:
        values = numeric(row.get(field) for row in rows)
        if values:
            result[field] = {
                "count": len(values), "mean": statistics.fmean(values),
                "median": statistics.median(values), "min": min(values),
                "max": max(values),
            }
    return result


def aggregate(rows: list[dict[str, Any]], keys: tuple[str, ...]) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[tuple(row.get(key) for key in keys)].append(row)
    return [{**dict(zip(keys, values)), **summarize(group)}
            for values, group in sorted(groups.items(), key=lambda item: str(item[0]))]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def selected_files(root: Path, include_predictions: bool) -> list[Path]:
    exact = (
        "analysis.json", "checksums.sha256", "scheduler.json", "status.json",
        "submitted-config.json", "job.log", "benchmark/manifest.json",
        "benchmark/summary.csv", "benchmark/summary.jsonl",
    )
    files = [root / name for name in exact if (root / name).is_file()]
    files.extend(path for path in (root / "system").glob("*") if path.is_file())
    artifact_names = {
        "result.json", "metrics.jsonl", "resource.json", "stderr.log",
        "stdout.log", "output.aln",
    }
    if include_predictions:
        artifact_names.add("prediction.fa")
    runs = root / "benchmark" / "runs"
    if runs.is_dir():
        files.extend(path for path in runs.rglob("*")
                     if path.is_file() and (
                         path.name in artifact_names or
                         (include_predictions and path.suffix in {".db", ".ct", ".bpp", ".pfs"})))
    return sorted(set(files), key=lambda path: str(path.relative_to(root)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--include-predictions", action="store_true")
    args = parser.parse_args()

    root = args.result_root.resolve()
    benchmark = root / "benchmark"
    rows = load_jsonl(benchmark / "summary.jsonl")
    result_files = sorted((benchmark / "runs").glob("*/*/rep-*/result.json"))
    analysis = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scheduler": load_json(root / "scheduler.json"),
        "job_status": load_json(root / "status.json"),
        "benchmark_manifest": load_json(benchmark / "manifest.json"),
        "integrity": {
            "summary_rows": len(rows), "result_json_files": len(result_files),
            "all_results_summarized": len(rows) == len(result_files),
            "unsuccessful_runs": [
                {key: row.get(key) for key in
                 ("condition", "dataset", "repetition", "status", "return_code")}
                for row in rows if row.get("status") != "success"
            ],
            "invalid_metric_runs": [
                {key: row.get(key) for key in ("condition", "dataset", "repetition")}
                for row in rows if row.get("metrics_valid") is not True
            ],
        },
        "overall": summarize(rows),
        "by_condition": aggregate(rows, ("condition",)),
        "by_condition_and_dataset": aggregate(rows, ("condition", "dataset")),
        "rows": rows,
    }
    analysis_path = root / "analysis.json"
    analysis_path.write_text(json.dumps(analysis, indent=2, sort_keys=True,
                                        allow_nan=False) + "\n", encoding="utf-8")

    files = selected_files(root, args.include_predictions)
    checksums_path = root / "checksums.sha256"
    checksum_lines = [f"{sha256(path)}  {path.relative_to(root)}" for path in files
                      if path != checksums_path]
    checksums_path.write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    files = selected_files(root, args.include_predictions)

    output = (args.output or (root / "analysis-bundle.tar.gz")).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(output, "w:gz") as archive:
        for path in files:
            if path.resolve() != output:
                archive.add(path, arcname=path.relative_to(root), recursive=False)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
