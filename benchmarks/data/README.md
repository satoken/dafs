# Benchmark datasets

Bulk inputs and references are generated locally and ignored by Git; the
versioned manifests record their provenance and selection rules.

```sh
python3 benchmarks/data/import_original.py
python3 benchmarks/data/prepare_latest_rfam.py
```

`import_original.py` reads the archives supplied in `~/dafs-data` and creates
`original-paper/manifest.json`: 691 PKfree, 82 PK, 52 MASTR, 1693 MXSCARNA,
and 18,990 BRAliBase 2.1 cases (21,508 runnable cases total). The BRAliBase
2.0 archive contains CT/structure annotations but not runnable sequence/MSA
pairs; it is recorded but cannot reproduce the paper's 599-case BRAliBase 2.0
experiment. The MASTR archive contains inputs without reference MSAs, so its
cases support runtime/SCI but not reference-based accuracy metrics.

`rfam15.1-current/manifest.json` is prepared from the official Rfam 15.1 seed
archive. It contains 545 families with usable consensus-structure tracks,
with a deterministic ten-sequence sample per family. Two seed families with
duplicate `SS_cons` tracks are excluded and named in the manifest rather than
silently resolved. The separate `rfam15.1-long/manifest.json` contains 19
reference-backed cases with at least one sequence of 500 nt or longer; each
ten-sequence sample is forced to include that family's longest sequence.
Eight cases include a sequence of at least 1,000 nt, and the longest is
4,065 nt. This is a useful long-RNA pilot, but not a broad or independent
long-RNA benchmark by itself.

`linearturbofold-nonviral/manifest.json` contains the 65 input groups shipped
with LinearTurboFold: 25 16S, 5 23S, 20 RNase P, 10 SRP, and 5 telomerase
groups. The longest sequence is 2,926 nt. Since this distribution has no
reference MSAs/structures for these groups, the benchmark reports runtime,
memory, output validation, and SCI, not alignment accuracy. In particular,
23S is a useful long-RNA scaling panel, but it should not be presented as a
reference-based accuracy result. The complete local run has 65/65 successful
cases; results are in `../results-linearturbofold-released-data/summary.csv`.

This 65-group repository collection is distinct from the LinearTurboFold
paper's accuracy experiment: the paper reports 60 five-sequence groups drawn
from four RNAStrAlign families and five resamplings. The supplied repository
does not include that exact group list or reference MSAs. The locally retained
RNAStrAlign archive contains per-sequence BPSEQ structures but no aligned
reference files, so it does not by itself reproduce the paper's alignment
accuracy evaluation.

The DAFS paper's original Rfam 11.0 selection list was not included in the
available files. `rfam11-paper-reconstructed/manifest.json` therefore remains
a clearly labeled reconstruction, not an exact reproduction. The supplied
local original-paper archives above are the exact files available from the
user, imported without replacing them with a reconstruction.

`murlet/manifest.json` and `rfam11-long/manifest.json` are separate established
and long-reference sets from the existing benchmark inventory. These remain
useful for complementary comparisons but are not substitutes for the
user-supplied DAFS paper archives.
