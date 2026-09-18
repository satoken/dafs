# TSUBAME CPU benchmark job

Build DAFS on a login node, prepare a benchmark JSON configuration, and submit
the job from the repository root. The script requests one `cpu_4` resource for
24 hours. Its optional final argument controls how many independent DAFS cases
run concurrently; use at most 4 with the current resource request.

```sh
qsub -g YOUR_TSUBAME_GROUP \
  benchmarks/tsubame/qsub_cpu.sh \
  benchmarks/configs/YOUR_CONFIG.json \
  /gs/bs/YOUR_GROUP/YOUR_NAME/dafs-results/RUN_NAME \
  build/src/dafs \
  benchmarks/tsubame/modules.example.sh \
  4
```

The binary, module-file, and worker-count arguments are optional; a direct
`qsub_cpu.sh` submission defaults to one worker. Use an empty string for the
module-file argument when specifying a worker count without a module file. Do not put
data preparation or compilation in this job: `module purge` is performed before
the optional module file is sourced, and the exact loaded modules are recorded.

The result directory is resumable. Re-submit the same command after a wall-time
or transient failure; completed cases are skipped by `benchmarks/run.py`.

The important outputs are:

- `analysis.json`: raw summary rows plus aggregates by condition and dataset.
- `analysis-bundle.tar.gz`: compact upload artifact for Codex analysis.
- `benchmark/summary.csv`: table for R, pandas, or spreadsheet analysis.
- `benchmark/summary.jsonl`: lossless line-oriented equivalent of the CSV.
- `benchmark/runs/.../metrics.jsonl`: per-iteration LB/UB and CBP metrics.
- `scheduler.json`, `status.json`, and `system/`: execution provenance.

Predicted alignments are retained in the result directory but excluded from the
default bundle to keep it small. To include them when diagnosing accuracy:

```sh
python3 benchmarks/tsubame/package_results.py RESULT_ROOT --include-predictions
```

## Linearization x RIBOSUM factorial benchmark

Submit the complete 2 x 2 comparison as one three-task array job:

```sh
qsub -g YOUR_TSUBAME_GROUP \
  benchmarks/tsubame/qsub_factorial.sh \
  /gs/bs/YOUR_GROUP/YOUR_NAME/dafs-results/factorial-RUN_NAME \
  build/src/dafs \
  benchmarks/tsubame/modules.example.sh
```

Task 1 runs 13 length-stratified Murlet dataset1 cases with references and
three repetitions. Task 2 runs all 25 nonviral LinearTurboFold 16S groups, and
task 3 runs all five 23S groups. Every suite compares the same four conditions:
CONTRAlign+CONTRAfold or LinearAlign+LinearPartition-C, each with RIBOSUM weight
0 or 0.075, with RNAalifold disabled. Each array task runs four cases
concurrently by default, so the array
can execute up to 12 single-threaded DAFS processes across three `cpu_4`
resources. Override the per-task worker count with a fourth argument after the
optional module file. The 16S and 23S nonlinear runs have finite per-case
timeouts so one quadratic/cubic case cannot consume the entire 24-hour
allocation.

## Nonlinear RIBOSUM weight sweep without RNAalifold

The dedicated Murlet sweep fixes the model to CONTRAlign + CONTRAfold with
`--no-alifold` and evaluates 15 RIBOSUM weights from 0 through 2.0. It runs the
13 length-stratified Murlet dataset1 cases once per weight (195 runs total) and
uses four workers by default:

```sh
qsub -g YOUR_TSUBAME_GROUP \
  benchmarks/tsubame/qsub_ribosum_sweep.sh \
  /gs/bs/YOUR_GROUP/YOUR_NAME/dafs-results/ribosum-sweep-001 \
  build/src/dafs \
  benchmarks/tsubame/modules.example.sh
```

Omit the module-file argument when no runtime modules are needed. The resulting
`analysis.json` and `benchmark/summary.csv` contain SPS, MCC, CBP F1, structure
sensitivity/PPV, runtime, memory, and internal optimization metrics for every
weight and dataset.
