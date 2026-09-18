#!/usr/bin/env bash
# Submit the Murlet RIBOSUM weight sweep as one parallel Slurm job.
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/../.." && pwd)

usage() {
    echo "Usage: $0 [RESULT_ROOT] [DAFS_BINARY] [VIENNA_VERSION] [JOBS] [NODELIST]" >&2
    echo "Runs 15 RIBOSUM weights for both nonlinear and linear DAFS on 13 Murlet cases." >&2
    echo "RNAalifold, RNAfold, and GNU time must already be available on PATH." >&2
    echo "NODELIST is optional and is passed to Slurm as --nodelist." >&2
}

if [[ $# -gt 5 ]]; then
    usage
    exit 2
fi

result_root=${1:-$repo_root/benchmarks/results-ribosum-sweep-murlet}
dafs_binary=${2:-$repo_root/build-amd64/dafs/src/dafs}
vienna_version=${3:-${VIENNA_VERSION:-2.5.1}}
jobs=${4:-4}
target_nodelist=${5:-}
config=$repo_root/benchmarks/configs/ribosum-sweep-murlet.json

absolute_path() {
    if [[ $1 = /* ]]; then
        printf '%s\n' "$1"
    else
        printf '%s/%s\n' "$repo_root" "$1"
    fi
}

result_root=$(absolute_path "$result_root")
dafs_binary=$(absolute_path "$dafs_binary")
mkdir -p "$result_root"

if [[ ! -f $config ]]; then
    echo "Configuration not found: $config" >&2
    exit 2
fi
if [[ ! -x $dafs_binary ]]; then
    echo "DAFS binary not executable: $dafs_binary" >&2
    exit 2
fi
if [[ ! $jobs =~ ^[1-9][0-9]*$ ]]; then
    echo "JOBS must be positive: $jobs" >&2
    exit 2
fi

sbatch_args=(
    --parsable
    --job-name=dafs-ribosum-sweep
    --export=ALL
    --cpus-per-task="$jobs"
    --mem=16G
    --time=24:00:00
    --output="$result_root/slurm-%x-%j.out"
    --error="$result_root/slurm-%x-%j.err"
)
if [[ -n $target_nodelist ]]; then
    sbatch_args+=(--nodelist="$target_nodelist")
fi

job_id=$(sbatch "${sbatch_args[@]}" \
    "$script_dir/run_ribosum_sweep.sbatch" \
    "$config" "$result_root" "$dafs_binary" "$vienna_version" "$jobs")

echo "dafs/ribosum-sweep-murlet: $job_id"
echo "Weights: 0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 2"
echo "Modes: nonlinear (CONTRAlign/CONTRAfold), linear (LinearAlign/lpc)"
echo "Results: $result_root"
if [[ -n $target_nodelist ]]; then
    echo "Slurm nodelist: $target_nodelist"
fi
