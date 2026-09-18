#!/usr/bin/env bash
#SBATCH --job-name=dafs-ribosum-sweep
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=dafs-ribosum-sweep-%j.out
#SBATCH --error=dafs-ribosum-sweep-%j.err

set -euo pipefail

usage() {
    echo "Usage: sbatch $0 CONFIG RESULT_ROOT DAFS_BINARY VIENNA_VERSION JOBS" >&2
}

if [[ $# -ne 5 ]]; then
    usage
    exit 2
fi

repo_root=${SLURM_SUBMIT_DIR:-$PWD}
exec bash "$repo_root/benchmarks/slurm/run_dafs.sbatch" "$@"
