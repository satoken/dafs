#!/bin/bash
#$ -cwd
#$ -l cpu_4=1
#$ -l h_rt=24:00:00
#$ -N dafs-benchmark
#$ -j y

set -euo pipefail
umask 027

usage() {
    echo "Usage: qsub -g GROUP $0 CONFIG RESULT_ROOT [BINARY] [MODULE_FILE] [JOBS]" >&2
    echo "Submit this script from the DAFS repository root." >&2
}

if [[ $# -lt 2 || $# -gt 5 ]]; then
    usage
    exit 2
fi

repo_root=${SGE_O_WORKDIR:-$PWD}
cd "$repo_root" || exit 2

absolute_path() {
    if [[ $1 = /* ]]; then
        printf '%s\n' "$1"
    else
        printf '%s/%s\n' "$repo_root" "$1"
    fi
}

config=$(absolute_path "$1")
result_root=$(absolute_path "$2")
binary=$(absolute_path "${3:-build/src/dafs}")
module_file=""
if [[ $# -ge 4 && -n $4 ]]; then
    module_file=$(absolute_path "$4")
fi
jobs=${5:-1}
if [[ ! $jobs =~ ^[1-9][0-9]*$ ]]; then
    echo "JOBS must be a positive integer (got $jobs)" >&2
    exit 2
fi
if (( jobs > 4 )); then
    echo "JOBS must not exceed the four CPUs requested by cpu_4=1 (got $jobs)" >&2
    exit 2
fi

if [[ ! -f $config ]]; then
    echo "Configuration not found: $config" >&2
    exit 2
fi
if [[ ! -x $binary ]]; then
    echo "DAFS binary not found or not executable: $binary" >&2
    exit 2
fi
if [[ -n $module_file && ! -f $module_file ]]; then
    echo "Module file not found: $module_file" >&2
    exit 2
fi

mkdir -p "$result_root/system"
result_root=$(cd "$result_root" && pwd)
job_id=${JOB_ID:-manual}
job_log="$result_root/job.log"
exec >"$job_log" 2>&1

export LC_ALL=C
export TZ=UTC
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
runner_status="not_started"

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
    echo "Bundle:   $result_root/analysis-bundle.tar.gz"
}
trap 'exit_code=$?; trap - EXIT; write_status_and_package "$exit_code"; exit "$exit_code"' EXIT

if type module >/dev/null 2>&1; then
    module purge
fi
if [[ -n $module_file ]]; then
    # The file should contain only module-load commands needed by this build.
    source "$module_file"
fi

cp -- "$config" "$result_root/submitted-config.json"

python3 - "$result_root/scheduler.json" "$config" "$binary" "$module_file" "$jobs" <<'PY'
import json
import os
import socket
import sys
from pathlib import Path

path, config, binary, module_file, jobs = sys.argv[1:]
keys = (
    "JOB_ID", "JOB_NAME", "HOSTNAME", "NSLOTS", "QUEUE", "SGE_O_HOST",
    "SGE_O_WORKDIR", "SGE_TASK_ID", "TMPDIR",
)
document = {
    "schema_version": 1,
    "scheduler": "Altair Grid Engine",
    "hostname": socket.gethostname(),
    "config": config,
    "binary": binary,
    "module_file": module_file or None,
    "benchmark_jobs": int(jobs),
    "environment": {key: os.environ.get(key) for key in keys},
    "thread_limits": {key: os.environ.get(key) for key in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    )},
}
Path(path).write_text(json.dumps(document, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
PY

uname -a >"$result_root/system/uname.txt" 2>&1 || true
lscpu >"$result_root/system/lscpu.txt" 2>&1 || true
cat /proc/meminfo >"$result_root/system/meminfo.txt" 2>&1 || true
if type module >/dev/null 2>&1; then
    module list >"$result_root/system/modules.txt" 2>&1 || true
fi
if command -v qstat >/dev/null 2>&1 && [[ $job_id != manual ]]; then
    qstat -j "$job_id" >"$result_root/system/qstat-job.txt" 2>&1 || true
fi
git rev-parse HEAD >"$result_root/system/git-commit.txt" 2>&1 || true
git status --short >"$result_root/system/git-status.txt" 2>&1 || true

echo "DAFS TSUBAME benchmark"
echo "Job:      $job_id"
echo "Host:     ${HOSTNAME:-unknown}"
echo "Config:   $config"
echo "Binary:   $binary"
echo "Results:  $result_root"
echo "Workers:  $jobs"
echo "Started:  $started_utc"

runner_status="running"
set +e
python3 -B "$repo_root/benchmarks/run.py" \
    --config "$config" \
    --binary "$binary" \
    --output-dir "$result_root/benchmark" \
    --jobs "$jobs"
runner_exit=$?
set -e
if [[ $runner_exit -eq 0 ]]; then
    runner_status="success"
else
    runner_status="failed"
fi
exit "$runner_exit"
