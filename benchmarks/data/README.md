# Benchmark datasets

Run `python3 benchmarks/data/prepare.py` from any directory to download and
verify the public source archives and generate DAFS inputs and references.
Downloaded and generated bulk files are ignored by Git; source metadata and
dataset manifests remain versionable.

`murlet/manifest.json` contains every alignment in the official Murlet paper
archive. Each aligned Stockholm reference is paired with an unaligned FASTA
input suitable for DAFS and can be loaded through `dataset_manifests` in a
benchmark configuration.

`raw/Rfam.seed.11.0.gz` is the official Rfam 11.0 seed used as the source for
the DAFS paper's 691 PKfree and 82 PK alignments. The original DAFS download
site is offline and its exact ten-sequence selections were not preserved in
the available archive. The seed is therefore retained as a verified raw source
and is not mislabeled as the exact DAFS PKfree/PK dataset.

`linearturbofold-nonviral/manifest.json` describes the nonviral sequence groups
distributed by the official LinearTurboFold repository. The 23S and 16S rRNA
groups are useful long-sequence scaling cases; RNase P, SRP, and telomerase are
included as shorter supporting cases. SARS-related collections are explicitly
excluded. The repository does not include the corresponding RNAStralign
reference alignments, so these inputs are for runtime/memory scaling and output
validation, not reference-based alignment accuracy scoring.
