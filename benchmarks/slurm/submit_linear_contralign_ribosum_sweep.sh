#!/usr/bin/env bash
# Submit LinearAlign-CONTRAlign RIBOSUM sweeps for Murlet, 16S, and 23S.
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/../.." && pwd)

usage() {
    echo "Usage: $0 [RESULT_ROOT] [DAFS_BINARY] [VIENNA_VERSION] [ARRAY_CONCURRENCY] [NODELIST] [VIENNA_BIN_DIR]" >&2
    echo "Runs the linear LinearAlign-CONTRAlign RIBOSUM sweep on Murlet, 16S, and 23S." >&2
    echo "Each array task processes one input FASTA and all 15 weights." >&2
    echo "RNAalifold, RNAfold, and GNU time must already be available on PATH." >&2
    echo "NODELIST is optional and is passed to Slurm as --nodelist." >&2
}

if [[ $# -gt 6 ]]; then
    usage
    exit 2
fi

result_root=${1:-$repo_root/benchmarks/results-linear-contralign-ribosum-sweep-murlet}
dafs_binary=${2:-$repo_root/build-amd64/dafs/src/dafs}
vienna_version=${3:-${VIENNA_VERSION:-2.5.1}}
array_concurrency=${4:-4}
target_nodelist=${5:-}
vienna_bin_dir=${6:-${VIENNA_BIN_DIR:-}}

absolute_path() {
    if [[ $1 = /* ]]; then
        printf '%s\n' "$1"
    else
        printf '%s/%s\n' "$repo_root" "$1"
    fi
}

result_root=$(absolute_path "$result_root")
dafs_binary=$(absolute_path "$dafs_binary")

if [[ -z $vienna_bin_dir && -n ${VIENNA_PREFIX:-} ]]; then
    vienna_bin_dir="$VIENNA_PREFIX/bin"
fi
if [[ -z $vienna_bin_dir ]] && command -v RNAfold >/dev/null 2>&1 &&
    command -v RNAalifold >/dev/null 2>&1; then
    vienna_bin_dir=$(dirname "$(command -v RNAfold)")
fi
if [[ -z $vienna_bin_dir ]]; then
    account_home=$(getent passwd "$(id -u)" 2>/dev/null | cut -d: -f6 || true)
    candidate="${account_home:-}/.local/app/viennarna-${vienna_version}/bin"
    if [[ -x $candidate/RNAfold && -x $candidate/RNAalifold ]]; then
        vienna_bin_dir=$candidate
    fi
fi
if [[ -n $vienna_bin_dir ]]; then
    vienna_bin_dir=$(absolute_path "$vienna_bin_dir")
fi

if [[ ! -x $dafs_binary ]]; then
    echo "DAFS binary not executable: $dafs_binary" >&2
    exit 2
fi
if [[ -z $vienna_bin_dir || ! -x $vienna_bin_dir/RNAfold || ! -x $vienna_bin_dir/RNAalifold ]]; then
    echo "ViennaRNA bin directory must contain executable RNAfold and RNAalifold." >&2
    echo "Set VIENNA_BIN_DIR or pass the sixth argument, for example:" >&2
    echo "  /home/sato-lab.org/satoken/.local/app/viennarna-${vienna_version}/bin" >&2
    exit 2
fi
if [[ ! $array_concurrency =~ ^[1-9][0-9]*$ ]]; then
    echo "ARRAY_CONCURRENCY must be positive: $array_concurrency" >&2
    exit 2
fi

dataset_count() {
    python3 - "$1" <<'PY'
import json
import sys
from pathlib import Path

config_path = Path(sys.argv[1]).resolve()
config = json.loads(config_path.read_text(encoding="utf-8"))
datasets = list(config.get("datasets", []))
for manifest_name in config.get("dataset_manifests", []):
    manifest_path = (config_path.parent / manifest_name).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    datasets.extend(manifest["datasets"])
dataset_filter = config.get("dataset_filter", {})
if dataset_filter.get("collections"):
    collections = {str(value) for value in dataset_filter["collections"]}
    datasets = [item for item in datasets
                if str(item.get("collection")) in collections]
if dataset_filter.get("ids"):
    ids = {str(value) for value in dataset_filter["ids"]}
    datasets = [item for item in datasets if str(item.get("id")) in ids]
if not datasets:
    raise SystemExit("configuration selected no datasets")
print(len(datasets))
PY
}

suites=(murlet 16s 23s)
configs=(
    "$repo_root/benchmarks/configs/linear-contralign-ribosum-sweep-murlet.json"
    "$repo_root/benchmarks/configs/linear-contralign-ribosum-sweep-16s.json"
    "$repo_root/benchmarks/configs/linear-contralign-ribosum-sweep-23s.json"
)
if [[ $result_root == *-murlet ]]; then
    output_roots=(
        "$result_root"
        "${result_root%-murlet}-16s"
        "${result_root%-murlet}-23s"
    )
else
    output_roots=(
        "$result_root/murlet"
        "$result_root/16s"
        "$result_root/23s"
    )
fi

for index in "${!suites[@]}"; do
    suite=${suites[$index]}
    config=${configs[$index]}
    output_dir=$(absolute_path "${output_roots[$index]}")

    if [[ ! -f $config ]]; then
        echo "Configuration not found: $config" >&2
        exit 2
    fi

    count=$(dataset_count "$config")
    concurrency=$array_concurrency
    if (( concurrency > count )); then
        concurrency=$count
    fi
    mkdir -p "$output_dir/tasks"

    array_args=(
        --parsable
        --job-name="dafs-linear-contralign-$suite"
        --array="0-$((count - 1))%$concurrency"
        --export=ALL
        --cpus-per-task=1
        --mem=16G
        --time=24:00:00
        --output="$output_dir/slurm-array-%A_%a.out"
        --error="$output_dir/slurm-array-%A_%a.err"
    )
    if [[ -n $target_nodelist ]]; then
        array_args+=(--nodelist="$target_nodelist")
    fi

    array_job_id=$(sbatch "${array_args[@]}" \
        "$script_dir/run_linear_contralign_ribosum_array.sbatch" \
        "$config" "$output_dir" "$dafs_binary" "$vienna_version" \
        "$vienna_bin_dir" "$repo_root")
    array_job_id=${array_job_id%%;*}

    merge_args=(
        --parsable
        --job-name="dafs-linear-contralign-$suite-merge"
        --dependency="afterany:$array_job_id"
        --export=ALL
        --cpus-per-task=1
        --mem=4G
        --time=01:00:00
        --output="$output_dir/slurm-merge-%j.out"
        --error="$output_dir/slurm-merge-%j.err"
    )
    if [[ -n $target_nodelist ]]; then
        merge_args+=(--nodelist="$target_nodelist")
    fi

    merge_job_id=$(sbatch "${merge_args[@]}" \
        "$script_dir/merge_linear_contralign_ribosum_array.sbatch" \
        "$config" "$output_dir" "$dafs_binary" "$vienna_version" \
        "$vienna_bin_dir" "$count" "$array_job_id" "$repo_root")

    echo "$suite: array $array_job_id ($count cases, max $concurrency concurrent), merge $merge_job_id"
    echo "  Results: $output_dir"
done

echo "LinearAlign-CONTRAlign RIBOSUM weights: 0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 2"
echo "ViennaRNA: $vienna_version"
echo "ViennaRNA bin: $vienna_bin_dir"
if [[ -n $target_nodelist ]]; then
    echo "Slurm nodelist: $target_nodelist"
fi
