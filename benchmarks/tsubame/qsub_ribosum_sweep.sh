#!/bin/bash
#$ -cwd
#$ -l cpu_4=1
#$ -l h_rt=24:00:00
#$ -N dafs-ribosum-sweep
#$ -j y

set -euo pipefail

if [[ $# -lt 1 || $# -gt 4 ]]; then
    echo "Usage: qsub -g GROUP $0 RESULT_ROOT [BINARY] [MODULE_FILE] [JOBS]" >&2
    exit 2
fi

repo_root=${SGE_O_WORKDIR:-$PWD}
result_root=$1
binary=${2:-build/src/dafs}
module_file=${3:-}
jobs=${4:-4}

exec bash "$repo_root/benchmarks/tsubame/qsub_cpu.sh" \
    "$repo_root/benchmarks/configs/ribosum-sweep-nonlinear-noalifold.json" \
    "$result_root" \
    "$binary" \
    "$module_file" \
    "$jobs"
