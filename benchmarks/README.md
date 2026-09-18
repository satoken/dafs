# DAFS benchmark infrastructure

The benchmark runner records predictions, GNU `time` resource measurements,
DAFS's internal JSONL trace, accuracy scores, an environment manifest, and
resumable per-run result files.  It uses only the Python standard library.

## Quick smoke run

```sh
python3 benchmarks/run.py \
  --config benchmarks/configs/smoke.json \
  --output-dir /tmp/dafs-benchmark-smoke
```

Completed cases are skipped when their command and input checksum match.
Use `--rerun`, `--rerun-failed`, or `--dry-run` to control execution.

For a 10-sequence example with five DD iterations per progressive merge:

```sh
python3 benchmarks/run.py --config benchmarks/configs/rf00005.json
```

## Configuration

The JSON configuration contains:

- `binary`: DAFS executable, relative to the configuration file;
- `output_dir`: result directory;
- `timeout_seconds`, `repetitions`, and `random_seed`;
- `environment`: variables pinned for every run;
- optional `sci`: a version-pinned ViennaRNA configuration. Use either
  `bin_dir` containing `RNAalifold` and `RNAfold`, or explicit `rnaalifold`
  and `rnafold` paths, together with the required `version` string and an
  optional per-tool `timeout_seconds`;
- `datasets`: objects with safe `id`, `input`, and optional `reference`;
- `dataset_manifests`: optional JSON files containing additional `datasets`;
- `dataset_filter`: optional `collections` and/or `ids` allowlists applied after
  loading datasets and manifests;
- `conditions`: objects with safe `id` and a list of DAFS `args`.

The maintained factorial configurations use the production defaults:
RIBOSUM weight `0.075` for enabled conditions and RNAalifold disabled. Use
`--alifold` only for an explicit ablation.

Dataset manifest paths are resolved relative to the manifest itself.  This
allows the eventual Murlet and PKfree archives to carry a versioned manifest
without putting hundreds of entries in an experiment configuration.

## Artifacts

Each case is stored below
`runs/<condition>/<dataset>/rep-XX/` and contains:

- `prediction.fa`: unchanged DAFS stdout;
- `stderr.log`: diagnostics;
- `metrics.jsonl`: stage and DD iteration events from DAFS;
- `resource.json`: wall/user/system time and maximum RSS from GNU `time`;
- `result.json`: checksums, command, resource data, invariant validation and
  accuracy metrics.

The runner regenerates `summary.jsonl` and `summary.csv` from the atomic
per-run results.  `manifest.json` records the binary checksum, Git revision,
configuration checksum, platform, Python version and controlled environment.

## Accuracy scorer

```sh
python3 benchmarks/score.py \
  --prediction prediction.fa \
  --reference reference.sto
```

References may be aligned FASTA with an `SS_cons` record or Stockholm with
`#=GC SS_cons`.  The scorer reports SPS, alignment PPV, structure sensitivity,
PPV, MCC, and pair-pair CBP precision/recall/F1.  Rfam/WUSS bracket symbols
are supported.  SCI is calculated when an explicit ViennaRNA installation is
supplied:

```sh
python3 benchmarks/score.py \
  --prediction prediction.fa \
  --vienna-bin-dir /opt/ViennaRNA-2.4.18/bin \
  --vienna-version 2.4.18
```

The SCI result records the consensus MFE, mean individual-sequence MFE, each
individual MFE, executable paths, and the reported ViennaRNA versions.  The
standard definition is `RNAalifold consensus MFE / mean RNAfold MFE`.

To enable SCI in `benchmarks/run.py`, add this to the benchmark JSON:

```json
"sci": {
  "bin_dir": "/opt/ViennaRNA-2.4.18/bin",
  "version": "2.4.18",
  "timeout_seconds": 900
}
```

For a site-dependent ViennaRNA location, the same configuration can instead
be overridden at launch with `--sci-bin-dir DIR --sci-version VERSION`.

The version is mandatory in benchmark configurations so results cannot be
silently mixed across ViennaRNA energy-model implementations. SCI can also be
computed for datasets without a reference alignment; in that case the result
contains SCI but the reference-based metrics remain absent.

## LinearTurboFold and Slurm comparison

`run_linearturbofold.py` evaluates the official LinearTurboFold output in the
same format as DAFS. It combines `output.aln` with the per-sequence `.db`
files, so structure sensitivity, PPV, MCC, CBP-F1, SPS, and SCI are available
when a reference is present. The nonviral LinearTurboFold manifests do not
contain reference structures; those cases therefore report runtime, memory,
output validation, and SCI only.

LinearTurboFold does not implement the DAFS RIBOSUM objective. It is run once
as `ribosum: not_applicable`; duplicating it as artificial RIBOSUM on/off
conditions would not be a meaningful ablation.

The Slurm entry point submits six independent jobs: DAFS and LinearTurboFold
for Murlet accuracy, 16S scaling, and 23S scaling. The DAFS jobs use the four
conditions nonlinear/linear × RIBOSUM off/on. The jobs use the same dataset
selection as the existing factorial configurations:

### Manual AMD64 build

Build both programs on an AMD64 node. Do not copy the ARM64 executables to the
cluster. Make GCC, CMake, pkg-config, ViennaRNA development files, and GLPK
available in your environment. The build needs the ViennaRNA installation
prefix, while the later SLURM jobs use the executables from `PATH`:

```sh
uname -m                         # must print x86_64/amd64
JOBS=8
VIENNA_PREFIX=/path/to/ViennaRNA
LINEARTURBOFOLD_SOURCE=/path/to/LinearTurboFold

export PATH="$VIENNA_PREFIX/bin:$PATH"
export PKG_CONFIG_PATH="$VIENNA_PREFIX/lib/pkgconfig:$VIENNA_PREFIX/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="$VIENNA_PREFIX/lib:$VIENNA_PREFIX/lib64:${LD_LIBRARY_PATH:-}"
pkg-config --modversion RNAlib2
command -v RNAalifold RNAfold

cmake -S . -B build-amd64/dafs \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/build-amd64/dafs-install" \
  -DCMAKE_PREFIX_PATH="$VIENNA_PREFIX" \
  -DBUILD_TESTING=OFF
cmake --build build-amd64/dafs --parallel "$JOBS"
cmake --install build-amd64/dafs

# -B forces recompilation, replacing any old ARM64 LinearTurboFold binary.
make -B -C "$LINEARTURBOFOLD_SOURCE" linearturbofold -j "$JOBS"

file build-amd64/dafs/src/dafs
file "$LINEARTURBOFOLD_SOURCE/bin/linearturbofold"
```

The expected outputs are `build-amd64/dafs/src/dafs` and
`$LINEARTURBOFOLD_SOURCE/bin/linearturbofold`. The DAFS CMake build uses the
GLPK backend by default when GLPK is available. `pkg-config --modversion
RNAlib2` must succeed.

After the manual build, submit the comparison jobs:

```sh
bash benchmarks/slurm/submit_comparison.sh \
  benchmarks/results-comparison \
  build-amd64/dafs/src/dafs \
  "$LINEARTURBOFOLD_SOURCE/bin/linearturbofold" \
  2.4.18 \
  4 \
  1 \
  node001
```

The last two numeric arguments are the number of concurrent benchmark cases
per DAFS job and per LinearTurboFold job. They are also passed to
`--cpus-per-task`. The submitter passes `--export=ALL` to each `sbatch`, so the
exported `PATH` (and `LD_LIBRARY_PATH`, if needed by a locally installed shared
library) is inherited by the batch job. The optional final `node001` argument
is passed to Slurm as `--nodelist=node001`; omit it to let the scheduler select
the execution node. `RERUN_BENCHMARKS=1` forces completed cases to be rerun;
otherwise matching completed cases are resumed.

The batch scripts verify that `RNAalifold` and `RNAfold` are discoverable on
`PATH` before starting. GNU `time` must also be available as an external
`time` command; if it is installed elsewhere, export `TIME_BINARY=/absolute/path/to/time`
before submission. The build uses the repository's CMake configuration for
DAFS and the official LinearTurboFold Makefile, and records source and binary
checksums in each benchmark manifest.

The six result directories are placed below the selected result root. Each
contains `summary.csv`, `analysis.json`, and a compressed analysis bundle; the
per-case `result.json` also records the SCI ViennaRNA version and executable
hashes. LinearTurboFold output identifiers are matched to the reference after
restoring the `/` separators that its FASTA writer changes to `_`.

If the method outputs already exist, `run_linearturbofold.py --score-only`
recalculates reference metrics and SCI without executing LinearTurboFold.
The Slurm wrapper exposes the same mode with `SCORE_ONLY=1`; it is useful for
repairing or extending evaluation after a scoring-only issue.

## Murlet RIBOSUM weight sweep

`ribosum-sweep-murlet.json` evaluates the 13 Murlet dataset1 cases used by the
accuracy benchmark at 15 weights (`0` through `2.0`) for both nonlinear
(CONTRAlign/CONTRAfold) and linear (LinearAlign/lpc) DAFS. RNAalifold is
disabled, and SCI is calculated with the version-pinned ViennaRNA tools on
`PATH`, for 390 runs in total:

```sh
bash benchmarks/slurm/submit_ribosum_sweep.sh \
  benchmarks/results-ribosum-sweep-murlet \
  build-amd64/dafs/src/dafs \
  2.5.1 \
  4 \
  node001
```

The final node argument is optional. The submitter passes the current `PATH`
to Slurm with `--export=ALL`; `RERUN_BENCHMARKS=1` reruns completed cases.
Results are grouped by the `nonlinear_ribosum_*` and `linear_ribosum_*`
conditions in `analysis.json`.

To rerun only the linear DAFS condition at RIBOSUM weight `0.3` for the
nonviral 16S and 23S scaling datasets (Murlet is already covered by the
weight sweep), submit both jobs with:

```sh
bash benchmarks/slurm/submit_linear_ribosum_03.sh \
  benchmarks/results-linear-ribosum-03 \
  build-amd64/dafs/src/dafs \
  2.5.1 \
  4 \
  node001
```

The 16S collection has 25 cases and the 23S collection has 5 cases. These
manifests do not contain reference structures, so these runs provide SCI,
runtime, memory, and output-validation results, but not SPS/MCC/reference-
based accuracy scores.

## LinearAlign-CONTRAlign RIBOSUM weight sweep

The submitter evaluates Murlet, 16S, and 23S with 15 weights using
`-a LinearAlign-CONTRAlign` and `lpc` folding. It creates one Slurm array per
collection, with one task per input FASTA; each task runs all 15 weights with
one worker. A dependent merge job collects each collection's task outputs into
`benchmark/summary.csv` and `analysis.json`:

```sh
bash benchmarks/slurm/submit_linear_contralign_ribosum_sweep.sh \
  benchmarks/results-linear-contralign-ribosum-sweep-murlet \
  build-amd64/dafs/src/dafs \
  2.5.1 \
  4 \
  node001 \
  /home/sato-lab.org/satoken/.local/app/viennarna-2.5.1/bin
```

The fifth argument is an optional Slurm nodelist and the sixth is the
ViennaRNA `bin` directory. The fourth argument is the maximum number of
simultaneous array tasks, not the number of workers inside a task. The
submitter passes the submitting environment with `--export=ALL` and the
explicit ViennaRNA directory to every array task, so `RNAalifold`, `RNAfold`,
and GNU `time` must be available on the selected compute nodes. The bin
directory is also auto-detected from `VIENNA_PREFIX`, `RNAfold`/`RNAalifold`
on the submission host, or the standard local-app path. Set `RERUN_BENCHMARKS=1`
before submission to rerun completed cases. Results are written below the given root as
`murlet/`, `16s/`, and `23s/`. For a root whose name ends in `-murlet`, the
Murlet results use that root itself and the 16S/23S results use sibling roots
ending in `-16s` and `-23s`, respectively.

## Internal trace

Pass `--metrics-jsonl FILE` directly to DAFS to collect:

- input size, models, thresholds, beams and sparse/dynamic modes;
- time and nnz for probability calculation and other major stages;
- per progressive merge and DD iteration: certified Lagrangian, best UB,
  feasible LB, gap, violations, Polyak scale and actual update, CBP counts,
  x/y/z beam scores, first-pruned-state and additive beam certificates,
  pruning counts, and
  intersection/consensus-repair feasible scores;
- aggregate DD time/iterations and CBP peak/add/pricing/removal counts.

The runner verifies `LB <= BestUB`, non-increasing BestUB, and non-decreasing
LB independently for every progressive merge.
