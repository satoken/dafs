#!/bin/bash
#$ -cwd
#$ -l cpu_4=1
#$ -l h_rt=24:00:00
#$ -N dafs-factorial
#$ -j y
#$ -t 1-3

set -euo pipefail

if [[ $# -lt 1 || $# -gt 4 ]]; then
    echo "Usage: qsub -g GROUP $0 RESULT_ROOT [BINARY] [MODULE_FILE] [JOBS]" >&2
    exit 2
fi

repo_root=${SGE_O_WORKDIR:-$PWD}
task_id=${SGE_TASK_ID:-0}
case "$task_id" in
    1) suite=factorial-accuracy ;;
    2) suite=factorial-16s ;;
    3) suite=factorial-23s ;;
    *) echo "SGE_TASK_ID must be 1, 2, or 3 (got $task_id)" >&2; exit 2 ;;
esac

result_root=$1
binary=${2:-build/src/dafs}
args=("$repo_root/benchmarks/configs/$suite.json" "$result_root/$suite" "$binary")
if [[ $# -ge 3 ]]; then
    args+=("$3")
else
    args+=("")
fi
jobs=${4:-4}
args+=("$jobs")

exec bash "$repo_root/benchmarks/tsubame/qsub_cpu.sh" "${args[@]}"
