#!/usr/bin/env bash
#SBATCH --job-name=dafs-linear-contralign-sweep
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=dafs-linear-contralign-sweep-%A_%a.out
#SBATCH --error=dafs-linear-contralign-sweep-%A_%a.err

set -euo pipefail
umask 027

usage() {
    echo "Usage: sbatch $0 CONFIG RESULT_ROOT DAFS_BINARY VIENNA_VERSION VIENNA_BIN_DIR REPO_ROOT" >&2
}

if [[ $# -ne 6 ]]; then
    usage
    exit 2
fi
if [[ -z ${SLURM_ARRAY_TASK_ID:-} ]]; then
    echo "This script must run as a Slurm job array" >&2
    exit 2
fi

vienna_bin_dir=$5
repo_root=$6
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
if [[ ! -x $vienna_bin_dir/RNAfold || ! -x $vienna_bin_dir/RNAalifold ]]; then
    echo "ViennaRNA executables not found in: $vienna_bin_dir" >&2
    exit 2
fi
export PATH="$vienna_bin_dir:${PATH:-}"
task_id=$SLURM_ARRAY_TASK_ID
task_root="$result_root/tasks/task-$task_id"
task_config="$task_root/config.json"

mkdir -p "$task_root"

dataset_id=$(python3 - "$config" "$task_config" "$task_id" <<'PY'
import json
import sys
from pathlib import Path

config_path = Path(sys.argv[1]).resolve()
task_config_path = Path(sys.argv[2]).resolve()
task_index = int(sys.argv[3])
config = json.loads(config_path.read_text(encoding="utf-8"))
config_base = config_path.parent

datasets = []
for dataset in config.get("datasets", []):
    resolved = dict(dataset)
    if "input" in resolved:
        input_path = Path(str(resolved["input"]))
        resolved["input"] = str((input_path if input_path.is_absolute()
                                  else config_base / input_path).resolve())
    if resolved.get("reference"):
        reference_path = Path(str(resolved["reference"]))
        resolved["reference"] = str((reference_path if reference_path.is_absolute()
                                      else config_base / reference_path).resolve())
    datasets.append(resolved)
for manifest_name in config.get("dataset_manifests", []):
    manifest_path = Path(str(manifest_name))
    if not manifest_path.is_absolute():
        manifest_path = config_base / manifest_path
    manifest_path = manifest_path.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_base = manifest_path.parent
    for dataset in manifest["datasets"]:
        resolved = dict(dataset)
        input_path = Path(str(resolved["input"]))
        resolved["input"] = str((input_path if input_path.is_absolute()
                                  else manifest_base / input_path).resolve())
        if resolved.get("reference"):
            reference_path = Path(str(resolved["reference"]))
            resolved["reference"] = str((reference_path if reference_path.is_absolute()
                                          else manifest_base / reference_path).resolve())
        datasets.append(resolved)

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
if task_index < 0 or task_index >= len(datasets):
    raise SystemExit(f"array index {task_index} is outside 0..{len(datasets)-1}")

selected = datasets[task_index]
selected_id = str(selected["id"])
config["datasets"] = [selected]
config.pop("dataset_manifests", None)
config.pop("dataset_filter", None)
config["output_dir"] = str((task_config_path.parent / "benchmark").resolve())
task_config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
print(selected_id)
PY
)

echo "Array task: $task_id"
echo "Dataset:    $dataset_id"
echo "Task root:  $task_root"

SLURM_SUBMIT_DIR="$repo_root" exec bash "$repo_root/benchmarks/slurm/run_dafs.sbatch" \
    "$task_config" "$task_root" "$binary" "$vienna_version" 1
