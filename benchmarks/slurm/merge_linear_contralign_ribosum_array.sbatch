#!/usr/bin/env bash
#SBATCH --job-name=dafs-linear-contralign-merge
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=dafs-linear-contralign-merge-%j.out
#SBATCH --error=dafs-linear-contralign-merge-%j.err

set -euo pipefail
umask 027

usage() {
    echo "Usage: sbatch $0 CONFIG RESULT_ROOT DAFS_BINARY VIENNA_VERSION VIENNA_BIN_DIR DATASET_COUNT ARRAY_JOB_ID REPO_ROOT" >&2
}

if [[ $# -ne 8 ]]; then
    usage
    exit 2
fi

vienna_bin_dir=$5
dataset_count=$6
array_job_id=$7
repo_root=$8
absolute_path() {
    if [[ $1 = /* ]]; then
        printf '%s\n' "$1"
    else
        printf '%s/%s\n' "$repo_root" "$1"
    fi
}

config=$(absolute_path "$1")
result_root=$(absolute_path "$2")
binary=$(absolute_path "$3")
vienna_version=$4

if [[ ! -d $repo_root ]]; then
    echo "Repository root not found: $repo_root" >&2
    exit 2
fi
if [[ ! $dataset_count =~ ^[1-9][0-9]*$ ]]; then
    echo "DATASET_COUNT must be positive: $dataset_count" >&2
    exit 2
fi

mkdir -p "$result_root"
exec >"$result_root/job.log" 2>&1

job_id=${SLURM_JOB_ID:-manual}
started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
runner_status=not_started

write_status() {
    local exit_code=$1
    local ended_utc
    ended_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    python3 - "$result_root/status.json" "$job_id" "$runner_status" \
        "$exit_code" "$started_utc" "$ended_utc" "$array_job_id" <<'PY'
import json
import sys
from pathlib import Path

path, job_id, runner_status, exit_code, started, ended, array_job_id = sys.argv[1:]
validation_path = Path(path).with_name("merge-validation.json")
Path(path).write_text(json.dumps({
    "schema_version": 1,
    "scheduler": "Slurm",
    "job_id": job_id,
    "array_job_id": array_job_id,
    "runner_status": runner_status,
    "exit_code": int(exit_code),
    "started_utc": started,
    "ended_utc": ended,
    "validation": (json.loads(validation_path.read_text(encoding="utf-8"))
                    if validation_path.is_file() else None),
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
}
trap 'exit_code=$?; if [[ $exit_code -ne 0 ]]; then runner_status=failed; fi; trap - EXIT; write_status "$exit_code"; exit "$exit_code"' EXIT

export LC_ALL=C
export TZ=UTC

runner_status=running
python3 - "$config" "$result_root" "$binary" "$vienna_version" "$vienna_bin_dir" "$dataset_count" "$array_job_id" <<'PY'
import csv
import json
import os
import shutil
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

config_path = Path(sys.argv[1]).resolve()
result_root = Path(sys.argv[2]).resolve()
binary = Path(sys.argv[3]).resolve()
vienna_version = sys.argv[4]
vienna_bin_dir = Path(sys.argv[5]).resolve()
dataset_count = int(sys.argv[6])
array_job_id = sys.argv[7]

task_roots = [result_root / "tasks" / f"task-{index}"
              for index in range(dataset_count)]
missing_summaries = []
for task_root in task_roots:
    if not (task_root / "benchmark" / "summary.jsonl").is_file():
        missing_summaries.append(str(task_root))
if missing_summaries:
    (result_root / "merge-validation.json").write_text(
        json.dumps({
            "array_job_id": array_job_id,
            "task_count": dataset_count,
            "missing_summaries": missing_summaries,
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    raise SystemExit("one or more array task summaries are missing")

merged_benchmark = result_root / "benchmark"
if merged_benchmark.exists():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = result_root / f"benchmark.previous.{stamp}"
    suffix = 1
    while backup.exists():
        backup = result_root / f"benchmark.previous.{stamp}.{suffix}"
        suffix += 1
    merged_benchmark.rename(backup)
merged_runs = merged_benchmark / "runs"
merged_runs.mkdir(parents=True, exist_ok=True)

rows = []
for task_root in task_roots:
    source_runs = task_root / "benchmark" / "runs"
    for source in sorted(source_runs.rglob("*")):
        if source.is_file():
            destination = merged_runs / source.relative_to(source_runs)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    summary_path = task_root / "benchmark" / "summary.jsonl"
    for line in summary_path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rows.append(json.loads(line))

rows.sort(key=lambda row: (str(row.get("condition")),
                           str(row.get("dataset")),
                           int(row.get("repetition", 0))))
(merged_benchmark / "summary.jsonl").write_text(
    "".join(json.dumps(row, sort_keys=True, allow_nan=False) + "\n"
            for row in rows), encoding="utf-8")
fieldnames = list(rows[0]) if rows else []
with (merged_benchmark / "summary.csv").open("w", encoding="utf-8", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

task_statuses = []
for task_root in task_roots:
    status_path = task_root / "status.json"
    task_statuses.append(
        json.loads(status_path.read_text(encoding="utf-8"))
        if status_path.is_file() else None)
validation = {
    "array_job_id": array_job_id,
    "task_count": dataset_count,
    "row_count": len(rows),
    "task_runner_failures": [
        {
            "task": index,
            "runner_status": (status or {}).get("runner_status"),
            "exit_code": (status or {}).get("exit_code"),
        }
        for index, status in enumerate(task_statuses)
        if status is None or status.get("runner_status") != "success"
    ],
    "unsuccessful_runs": [
        {key: row.get(key) for key in
         ("condition", "dataset", "repetition", "status", "return_code")}
        for row in rows if row.get("status") != "success"
    ],
    "invalid_metric_runs": [
        {key: row.get(key) for key in ("condition", "dataset", "repetition")}
        for row in rows if row.get("metrics_valid") is not True
    ],
}
(result_root / "merge-validation.json").write_text(
    json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")

first_manifest = json.loads(
    (task_roots[0] / "benchmark" / "manifest.json").read_text(encoding="utf-8"))
first_manifest["array_job_id"] = array_job_id
first_manifest["array_task_count"] = dataset_count
first_manifest["config_path"] = str(config_path)
first_manifest["binary"] = str(binary)
first_manifest["sci_version_override"] = vienna_version
first_manifest["sci_bin_dir_override"] = str(vienna_bin_dir)
(merged_benchmark / "manifest.json").write_text(
    json.dumps(first_manifest, indent=2, sort_keys=True) + "\n",
    encoding="utf-8")

result_root.joinpath("submitted-config.json").write_text(
    config_path.read_text(encoding="utf-8"), encoding="utf-8")
keys = ("SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_SUBMIT_DIR",
        "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_COUNT",
        "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE",
        "SLURM_JOB_PARTITION")
(result_root / "scheduler.json").write_text(
    json.dumps({
        "schema_version": 1,
        "scheduler": "Slurm",
        "hostname": socket.gethostname(),
        "config": str(config_path),
        "binary": str(binary),
        "vienna_version": vienna_version,
        "vienna_bin_dir": str(vienna_bin_dir),
        "benchmark_jobs": 1,
        "array_job_id": array_job_id,
        "array_task_count": dataset_count,
        "environment": {key: os.environ.get(key) for key in keys},
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
system = result_root / "system"
system.mkdir(exist_ok=True)
PY

uname -a >"$result_root/system/uname.txt" 2>&1 || true
lscpu >"$result_root/system/lscpu.txt" 2>&1 || true
git -C "$repo_root" rev-parse HEAD >"$result_root/system/git-commit.txt" 2>&1 || true
git -C "$repo_root" status --short >"$result_root/system/git-status.txt" 2>&1 || true

if python3 - "$result_root/merge-validation.json" <<'PY'
import json
import sys

validation = json.load(open(sys.argv[1], encoding="utf-8"))
failed = (validation["task_runner_failures"] or
          validation["unsuccessful_runs"] or
          validation["invalid_metric_runs"])
raise SystemExit(1 if failed else 0)
PY
then
    runner_status=success
    runner_exit=0
else
    runner_status=failed
    runner_exit=1
fi

# Write status before packaging so analysis.json records the merge result.
write_status "$runner_exit"
python3 "$repo_root/benchmarks/tsubame/package_results.py" "$result_root"
echo "Array job: $array_job_id"
echo "Summary:   $result_root/benchmark/summary.csv"
echo "Analysis:  $result_root/analysis.json"
exit "$runner_exit"
