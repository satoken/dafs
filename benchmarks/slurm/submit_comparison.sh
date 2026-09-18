#!/usr/bin/env bash
# Submit DAFS and LinearTurboFold suites independently so Slurm can run them
# concurrently. Each suite internally runs JOBS cases in parallel.
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/../.." && pwd)

usage() {
    echo "Usage: $0 [RESULT_ROOT] [DAFS_BINARY] [LINEARTURBOFOLD_BINARY] [VIENNA_VERSION] [DAFS_JOBS] [LTF_JOBS] [NODELIST]" >&2
    echo "RNAalifold and RNAfold must already be available on PATH." >&2
    echo "NODELIST is an optional Slurm node name, passed as --nodelist." >&2
}

if [[ $# -gt 7 ]]; then usage; exit 2; fi

result_root=${1:-$repo_root/benchmarks/results-comparison}
dafs_binary=${2:-$repo_root/build-amd64/dafs/src/dafs}
ltf_binary=${3:-$repo_root/../../LinearFold/LinearTurboFold/bin/linearturbofold}
vienna_version=${4:-${VIENNA_VERSION:-2.4.18}}
dafs_jobs=${5:-4}
ltf_jobs=${6:-1}
target_nodelist=${7:-}

absolute_path() {
    if [[ $1 = /* ]]; then
        printf '%s\n' "$1"
    else
        printf '%s/%s\n' "$repo_root" "$1"
    fi
}

result_root=$(absolute_path "$result_root")
dafs_binary=$(absolute_path "$dafs_binary")
ltf_binary=$(absolute_path "$ltf_binary")
mkdir -p "$result_root"

submit_one() {
    local method=$1
    local suite=$2
    local config=$3
    local binary=$4
    local jobs=$5
    local memory=$6
    local script=$7
    local job_id
    local common_args=("$config" "$result_root/$method-$suite" "$binary"
                       "$vienna_version" "$jobs")
    local sbatch_args=(
        --parsable
        --job-name="${method}-${suite}"
        --export=ALL
        --cpus-per-task="$jobs"
        --mem="$memory"
        --output="$result_root/slurm-%x-%j.out"
        --error="$result_root/slurm-%x-%j.err"
    )
    if [[ -n $target_nodelist ]]; then
        sbatch_args+=(--nodelist="$target_nodelist")
    fi
    job_id=$(sbatch "${sbatch_args[@]}" \
        "$script" "${common_args[@]}")
    echo "$method/$suite: $job_id"
}

echo "Submitting comparison suites under $result_root"
echo "DAFS RIBOSUM conditions: nonlinear/linear x on/off"
echo "LinearTurboFold RIBOSUM: not applicable (one baseline condition)"
if [[ -n $target_nodelist ]]; then echo "Slurm nodelist: $target_nodelist"; fi

submit_one dafs accuracy \
    "$repo_root/benchmarks/configs/factorial-accuracy.json" \
    "$dafs_binary" "$dafs_jobs" 16G "$script_dir/run_dafs.sbatch"
submit_one dafs 16s \
    "$repo_root/benchmarks/configs/factorial-16s.json" \
    "$dafs_binary" "$dafs_jobs" 16G "$script_dir/run_dafs.sbatch"
submit_one dafs 23s \
    "$repo_root/benchmarks/configs/factorial-23s.json" \
    "$dafs_binary" "$dafs_jobs" 16G "$script_dir/run_dafs.sbatch"

submit_one linearturbofold accuracy \
    "$repo_root/benchmarks/configs/linearturbofold-factorial-accuracy.json" \
    "$ltf_binary" "$ltf_jobs" 32G "$script_dir/run_linearturbofold.sbatch"
submit_one linearturbofold 16s \
    "$repo_root/benchmarks/configs/linearturbofold-factorial-16s.json" \
    "$ltf_binary" "$ltf_jobs" 32G "$script_dir/run_linearturbofold.sbatch"
submit_one linearturbofold 23s \
    "$repo_root/benchmarks/configs/linearturbofold-factorial-23s.json" \
    "$ltf_binary" "$ltf_jobs" 32G "$script_dir/run_linearturbofold.sbatch"
