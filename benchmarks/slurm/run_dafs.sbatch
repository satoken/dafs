#!/usr/bin/env bash
#SBATCH --job-name=dafs-benchmark
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --output=dafs-benchmark-%j.out
#SBATCH --error=dafs-benchmark-%j.err

set -euo pipefail
umask 027

usage() {
    echo "Usage: sbatch $0 CONFIG RESULT_ROOT DAFS_BINARY VIENNA_VERSION JOBS" >&2
}

if [[ $# -ne 5 ]]; then
    usage
    exit 2
fi

repo_root=${SLURM_SUBMIT_DIR:-$PWD}
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
jobs=$5

case "$(uname -m)" in
    x86_64|amd64) ;;
    *)
        echo "SLURM job is not running on AMD64: $(uname -m)" >&2
        exit 1
        ;;
esac
if [[ ! -f $config ]]; then echo "Configuration not found: $config" >&2; exit 2; fi
if [[ ! -x $binary ]]; then echo "DAFS binary not executable: $binary" >&2; exit 2; fi
rnaalifold=$(command -v RNAalifold || true)
rnafold=$(command -v RNAfold || true)
if [[ -z $rnaalifold || -z $rnafold ]]; then
    echo "RNAalifold and RNAfold must be available on PATH" >&2
    exit 2
fi
if [[ ! $jobs =~ ^[1-9][0-9]*$ ]]; then echo "JOBS must be positive" >&2; exit 2; fi
if (( jobs > ${SLURM_CPUS_PER_TASK:-1} )); then
    echo "JOBS=$jobs exceeds SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}" >&2
    exit 2
fi

mkdir -p "$result_root/system"
result_root=$(cd "$result_root" && pwd)
job_id=${SLURM_JOB_ID:-manual}
started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
runner_status=not_started
exec >"$result_root/job.log" 2>&1

write_status_and_package() {
    local exit_code=$1
    local ended_utc
    set +e
    ended_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    python3 - "$result_root/status.json" "$job_id" "$runner_status" \
        "$exit_code" "$started_utc" "$ended_utc" <<'PY'
import json
import sys
from pathlib import Path

path, job_id, runner_status, exit_code, started, ended = sys.argv[1:]
Path(path).write_text(json.dumps({
    "schema_version": 1,
    "scheduler": "Slurm",
    "job_id": job_id,
    "runner_status": runner_status,
    "exit_code": int(exit_code),
    "started_utc": started,
    "ended_utc": ended,
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
    python3 "$repo_root/benchmarks/tsubame/package_results.py" "$result_root" || true
    echo "Exit code: $exit_code"
    echo "Analysis: $result_root/analysis.json"
}
trap 'exit_code=$?; trap - EXIT; write_status_and_package "$exit_code"; exit "$exit_code"' EXIT

export LC_ALL=C
export TZ=UTC
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cp -- "$config" "$result_root/submitted-config.json"
uname -a >"$result_root/system/uname.txt" 2>&1 || true
lscpu >"$result_root/system/lscpu.txt" 2>&1 || true
cat /proc/meminfo >"$result_root/system/meminfo.txt" 2>&1 || true
if type module >/dev/null 2>&1; then module list >"$result_root/system/modules.txt" 2>&1 || true; fi
git -C "$repo_root" rev-parse HEAD >"$result_root/system/git-commit.txt" 2>&1 || true
git -C "$repo_root" status --short >"$result_root/system/git-status.txt" 2>&1 || true

python3 - "$result_root/scheduler.json" "$config" "$binary" "$rnaalifold" "$rnafold" \
    "$vienna_version" "$jobs" <<'PY'
import json
import os
import socket
import sys
from pathlib import Path

path, config, binary, rnaalifold, rnafold, vienna_version, jobs = sys.argv[1:]
keys = ("SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_SUBMIT_DIR",
        "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_JOB_PARTITION")
Path(path).write_text(json.dumps({
    "schema_version": 1,
    "scheduler": "Slurm",
    "hostname": socket.gethostname(),
    "config": config,
    "binary": binary,
    "rnaalifold": rnaalifold,
    "rnafold": rnafold,
    "vienna_version": vienna_version,
    "benchmark_jobs": int(jobs),
    "environment": {key: os.environ.get(key) for key in keys},
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

echo "DAFS benchmark"
echo "Job:      $job_id"
echo "Host:     ${HOSTNAME:-unknown} ($(uname -m))"
echo "Config:   $config"
echo "Binary:   $binary"
echo "Results:  $result_root"
echo "Workers:  $jobs"
echo "SCI:      ViennaRNA $vienna_version"

runner_status=running
set +e
rerun_args=()
if [[ ${RERUN_BENCHMARKS:-0} == 1 ]]; then rerun_args+=(--rerun); fi
python3 -B "$repo_root/benchmarks/run.py" \
    --config "$config" \
    --binary "$binary" \
    --output-dir "$result_root/benchmark" \
    --jobs "$jobs" \
    --sci-version "$vienna_version" \
    "${rerun_args[@]}"
runner_exit=$?
set -e
if [[ $runner_exit -eq 0 ]]; then runner_status=success; else runner_status=failed; fi
exit "$runner_exit"
