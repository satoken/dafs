#!/usr/bin/env bash
# Submit linear DAFS with RIBOSUM weight 0.3 for the 16S and 23S datasets.
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/../.." && pwd)

usage() {
    echo "Usage: $0 [RESULT_ROOT] [DAFS_BINARY] [VIENNA_VERSION] [JOBS] [NODELIST]" >&2
    echo "Submits separate 16S and 23S jobs for linear DAFS at RIBOSUM weight 0.3." >&2
    echo "RNAalifold, RNAfold, and GNU time must already be available on PATH." >&2
    echo "NODELIST is optional and is passed to Slurm as --nodelist." >&2
}

if [[ $# -gt 5 ]]; then
    usage
    exit 2
fi

result_root=${1:-$repo_root/benchmarks/results-linear-ribosum-03}
dafs_binary=${2:-$repo_root/build-amd64/dafs/src/dafs}
vienna_version=${3:-${VIENNA_VERSION:-2.5.1}}
jobs=${4:-4}
target_nodelist=${5:-}

absolute_path() {
    if [[ $1 = /* ]]; then
        printf '%s\n' "$1"
    else
        printf '%s/%s\n' "$repo_root" "$1"
    fi
}

result_root=$(absolute_path "$result_root")
dafs_binary=$(absolute_path "$dafs_binary")

if [[ ! -x $dafs_binary ]]; then
    echo "DAFS binary not executable: $dafs_binary" >&2
    exit 2
fi
if ! [[ $jobs =~ ^[1-9][0-9]*$ ]]; then
    echo "JOBS must be positive: $jobs" >&2
    exit 2
fi

for collection in 16s 23s; do
    config="$repo_root/benchmarks/configs/linear-ribosum-03-${collection}.json"
    if [[ ! -f $config ]]; then
        echo "Configuration not found: $config" >&2
        exit 2
    fi
done

submit_one() {
    local collection=$1
    local config="$repo_root/benchmarks/configs/linear-ribosum-03-${collection}.json"
    local output_dir="$result_root/${collection}"
    local job_id
    mkdir -p "$output_dir"

    local sbatch_args=(
        --parsable
        --job-name="dafs-linear-ribosum-03-${collection}"
        --export=ALL
        --cpus-per-task="$jobs"
        --mem=16G
        --time=24:00:00
        --output="$output_dir/slurm-%x-%j.out"
        --error="$output_dir/slurm-%x-%j.err"
    )
    if [[ -n $target_nodelist ]]; then
        sbatch_args+=(--nodelist="$target_nodelist")
    fi

    job_id=$(sbatch "${sbatch_args[@]}" \
        "$script_dir/run_dafs.sbatch" \
        "$config" "$output_dir" "$dafs_binary" "$vienna_version" "$jobs")
    echo "linear-ribosum-03/${collection}: $job_id"
}

echo "Submitting linear DAFS at RIBOSUM weight 0.3"
echo "Results: $result_root"
echo "Workers per job: $jobs"
echo "ViennaRNA: $vienna_version"
if [[ -n $target_nodelist ]]; then
    echo "Slurm nodelist: $target_nodelist"
fi

submit_one 16s
submit_one 23s
