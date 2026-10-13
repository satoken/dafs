# External structural-alignment tools

Pins, source checksums, and build/runtime status are recorded in
`versions.json`. The current comparison configs use LocARNA 2.0.1,
RNAmountAlign 1.0, LaRA2 2.0.1 with T-Coffee 13.46.2, RAF 1.00,
MAFFT 7.526 Q-INS-i/X-INS-i/G-INS-i, and the pinned LinearTurboFold source.
For RAF, the archived official source of CONTRAlign 2.01 is also built; RAF
uses its `--vienna` mode with the locally pinned ViennaRNA 2.7.2.

The source and install trees under `benchmarks/tools/` are ignored build
artifacts. These commands describe the source pins and legacy compiler
compatibility flags used here; use a Linux environment with GCC, Make, Perl,
and ViennaRNA 2.7.2 available.

## ViennaRNA 2.7.2

RNAmountAlign and LaRA2 are linked against the same local ViennaRNA build:

```sh
curl -fL 'https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.2/ViennaRNA-2.7.2.tar.gz' \
  -o benchmarks/tools/downloads/ViennaRNA-2.7.2.tar.gz
tar -xzf benchmarks/tools/downloads/ViennaRNA-2.7.2.tar.gz \
  -C benchmarks/tools/sources
VIENNA_PREFIX="$PWD/benchmarks/tools/install/viennarna-2.7.2"
cd benchmarks/tools/sources/ViennaRNA-2.7.2
./configure --prefix="$VIENNA_PREFIX" --without-swig --without-python \
  --without-kinfold --without-forester --without-rnalocmin --without-gsl \
  --disable-lto
make -j4
make install
```

## LocARNA 2.0.1

Source tag `v2.0.1`, commit `a6736722216bbb9541caf397f87742d5ff099d0b`:

```sh
git clone https://github.com/s-will/LocARNA benchmarks/tools/sources/LocARNA-2.0.1
git -C benchmarks/tools/sources/LocARNA-2.0.1 checkout --detach v2.0.1
cd benchmarks/tools/sources/LocARNA-2.0.1
autoreconf -i
./configure --prefix="$PWD/../../install/locarna-2.0.1" \
  --with-vrna=/path/to/ViennaRNA-2.7.2
make -C src clean
make -C src -j4 CXXFLAGS='-O3 -DNDEBUG -fsigned-char'
make -C src install CXXFLAGS='-O3 -DNDEBUG -fsigned-char'
make -C lib install
```

On this AArch64 host, the upstream option parser stores `getopt`'s EOF result
in a plain `char`; since GCC defaults plain char to unsigned here, LocARNA
loops during option parsing unless it is compiled with `-fsigned-char`. Build
and therefore include `-fsigned-char` in both build and install invocations.

The full documentation build needs `help2man`; executable targets do not.
LocARNA 2.0.1 passed the one-family Rfam smoke run after this compatibility
flag was applied.

## RAF 1.00 and its alignment dependency

The final RAF 1.00 source archive is recovered from the Internet Archive's
capture of the Stanford distribution. CONTRAlign 2.01 is required by RAF even
in `--vienna` folding mode.

```sh
mkdir -p benchmarks/tools/sources/raf benchmarks/tools/install/raf-1.00/bin
curl -fL 'https://web.archive.org/web/20100813092228id_/http://contra.stanford.edu/contrafold/raf_v1_00.tar.gz' \
  -o benchmarks/tools/downloads/raf_v1_00.tar.gz
tar -xzf benchmarks/tools/downloads/raf_v1_00.tar.gz -C benchmarks/tools/sources
cd benchmarks/tools/sources/raf/src
make -j4 CXX='g++ -include climits'
install -m 755 raf ../../../../install/raf-1.00/bin/raf

curl -fL 'https://web.archive.org/web/20100813141936id_/http://contra.stanford.edu/contralign/contralign_v2_01.tar.gz' \
  -o ../../downloads/contralign_v2_01.tar.gz
mkdir -p ../contralign-2.01
tar -xzf ../../downloads/contralign_v2_01.tar.gz -C ../contralign-2.01 --strip-components=1
cd ../contralign-2.01/src
make -j4 MODEL_TYPE=-DRNA=1 CXX='g++ -include climits -fpermissive'
install -Dm755 contralign ../../../../install/contralign-2.01/bin/contralign
```

RAF's comparison config points `VIENNA_RNA_DIR` at the ViennaRNA `bin` folder,
`CONTRALIGN_DIR` at the built CONTRAlign `bin` folder, and uses
`raf predict --vienna INPUT`.

## MAFFT 7.526

The official GitLab `v7.526` source includes the structural extensions needed
for Q-INS-i. Q-INS-i and X-INS-i are invoked through their method-specific
symlink names; the runner deliberately preserves those symlinks.

```sh
git clone https://gitlab.com/sysimm/mafft.git benchmarks/tools/sources/mafft-v7.526-gitlab
git -C benchmarks/tools/sources/mafft-v7.526-gitlab checkout --detach v7.526
MAFFT_PREFIX="$PWD/benchmarks/tools/install/mafft-7.526"
make -C benchmarks/tools/sources/mafft-v7.526-gitlab/core -j4
make -C benchmarks/tools/sources/mafft-v7.526-gitlab/core install PREFIX="$MAFFT_PREFIX"
make -C benchmarks/tools/sources/mafft-v7.526-gitlab/extensions -j4
make -C benchmarks/tools/sources/mafft-v7.526-gitlab/extensions install PREFIX="$MAFFT_PREFIX"
```

The configs set `MAFFT_BINARIES` relative to their own file location so
MAFFT's helper programs continue to resolve after the runner changes cwd.

## RNAmountAlign 1.0

RNAmountAlign is pinned to upstream commit
`2a9579ecac6cb612dfb12231eb99be502d736ee8`. The distribution bundles an
x86-only ViennaRNA static library; build against the local ViennaRNA source
instead, and apply the small linker fix for an upstream `verbose` symbol:

```sh
git clone https://github.com/clotelab/RNAmountAlign benchmarks/tools/sources/RNAmountAlign-1.0
git -C benchmarks/tools/sources/RNAmountAlign-1.0 checkout --detach 2a9579ecac6cb612dfb12231eb99be502d736ee8
git -C benchmarks/tools/sources/RNAmountAlign-1.0 apply "$PWD/benchmarks/tools/patches/rnamountalign-v1.0-linker.patch"
VIENNA_PREFIX="$PWD/benchmarks/tools/install/viennarna-2.7.2"
make -C benchmarks/tools/sources/RNAmountAlign-1.0/src -j4 \
  LDFLAGS="-L$VIENNA_PREFIX/lib -L. -lm -lRNA -fopenmp -pthread"
mkdir -p benchmarks/tools/install/rnamountalign-1.0/bin
install -m 755 benchmarks/tools/sources/RNAmountAlign-1.0/src/RNAmountAlign \
  benchmarks/tools/install/rnamountalign-1.0/bin/RNAmountAlign
cp -R benchmarks/tools/sources/RNAmountAlign-1.0/src/matrices \
  benchmarks/tools/sources/RNAmountAlign-1.0/src/energy_params \
  benchmarks/tools/install/rnamountalign-1.0/bin/
```

The benchmark invokes progressive multiple alignment in FASTA mode. Its input
alphabet is limited to A/C/G/U/T; inputs with other symbols are recorded as
`skipped_incompatible_input`. Symbols are never rewritten or removed.

## LaRA2 2.0.1 with T-Coffee 13.46.2

LaRA2 computes pairwise structural alignments and emits a T-Coffee library; it
does not itself produce a final MSA. The configured method follows the upstream
workflow: LaRA2 library output is passed to T-Coffee, which builds the MSA.
This is separate from the native MAFFT Q-INS-i/X-INS-i/G-INS-i baselines.

LaRA2 is pinned to upstream commit
`f94301fdcb61db59025b451299f8be9aff4fd94f` (SeqAn 2.4.1). Build it against
the local ViennaRNA installation:

```sh
git clone --recurse-submodules https://github.com/seqan/lara benchmarks/tools/sources/LaRA2-v2.0.1
git -C benchmarks/tools/sources/LaRA2-v2.0.1 checkout --detach f94301fdcb61db59025b451299f8be9aff4fd94f
git -C benchmarks/tools/sources/LaRA2-v2.0.1 submodule update --init --recursive
VIENNA_PREFIX="$PWD/benchmarks/tools/install/viennarna-2.7.2"
LARA_PREFIX="$PWD/benchmarks/tools/install/lara-2.0.1"
cmake -S benchmarks/tools/sources/LaRA2-v2.0.1 \
  -B benchmarks/tools/sources/LaRA2-v2.0.1/build \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_RUNTIME_OUTPUT_DIRECTORY="$LARA_PREFIX/bin" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DLARA_COMPILE_THREADS=4 \
  -DCMAKE_CXX_FLAGS="-I$VIENNA_PREFIX/include" \
  -DVIENNA_RNA_LIB="$VIENNA_PREFIX/lib/libRNA.a" \
  -DVIENNA_RNA_PATH="$VIENNA_PREFIX/include"
cmake --build benchmarks/tools/sources/LaRA2-v2.0.1/build -j4
```

T-Coffee is pinned to upstream commit
`04ec7406b2e13cd18c3f2077541079c8e059ff01` (version
`13.46.2.7c9e712d`):

```sh
git clone https://github.com/cbcrg/tcoffee.git benchmarks/tools/sources/tcoffee
git -C benchmarks/tools/sources/tcoffee checkout --detach 04ec7406b2e13cd18c3f2077541079c8e059ff01
make -C benchmarks/tools/sources/tcoffee/t_coffee/src t_coffee -j4 CFLAGS='-O2 -w'
mkdir -p benchmarks/tools/install/tcoffee/bin
install -m755 benchmarks/tools/sources/tcoffee/t_coffee/src/t_coffee \
  benchmarks/tools/install/tcoffee/bin/t_coffee
```

The common-runner wrapper is `lara2_tcoffee_wrapper.py`. It writes the LaRA2
library and final MSA in each case's `tool/` directory, and gives T-Coffee
case-local temporary and home directories so it does not write into the user's
home. Configs set `LARA2_BIN_DIR` and `TCOFFEE_BIN_DIR`.

## LinearTurboFold

The official repository is pinned to
`445808976857df01c57e5d9302a0bbee10515bb5`:

```sh
git clone https://github.com/LinearFold/LinearTurboFold benchmarks/tools/sources/LinearTurboFold
git -C benchmarks/tools/sources/LinearTurboFold checkout --detach 445808976857df01c57e5d9302a0bbee10515bb5
make -C benchmarks/tools/sources/LinearTurboFold linearturbofold -j4
```

Run the official 65-group nonviral dataset using its specialized parser and
scorer:

```sh
python3 benchmarks/run_linearturbofold.py \
  --config benchmarks/configs/linearturbofold-released-data.json \
  --sci-bin-dir /path/to/ViennaRNA-2.7.2/bin --sci-version 2.7.2 --jobs 2
```

The same groups can be run through the available structural-alignment tools
with `competitors-linearturbofold-data.json`. No accuracy score should be
reported for that input-only collection because the repository does not ship
reference MSAs or structures.

## CentroidAlign source availability

The DAFS paper's CentroidAlign baseline is the historical version 1.00. Its
official page still appears in archived records, but the linked source archive
is not present in the captures and the live download endpoint is unavailable.
The related CentroidAlifold package is not the CentroidAlign algorithm and is
not substituted. Consequently CentroidAlign is explicitly marked unavailable
in `versions.json` and is omitted from executable run configs until its actual
source archive is recovered.

## Current-tool benchmark commands

The Rfam 15.1, long-RNA, original-archive, and LinearTurboFold dataset configs
are under `benchmarks/configs/competitors-*.json`. Start with the small smoke
case, which checks RNAmountAlign, LaRA2+T-Coffee, RAF, the three MAFFT modes,
output normalization, accuracy, and SCI:

```sh
python3 benchmarks/run_competitor.py \
  --config benchmarks/configs/competitors-smoke-rfam15.1.json \
  --sci-bin-dir /path/to/ViennaRNA-2.7.2/bin --sci-version 2.7.2 --jobs 2
```

LocARNA, RNAmountAlign, LaRA2+T-Coffee, RAF, and the three MAFFT modes are wired
to the common runner. Full original-archive runs cover 21,508 cases and can be
long; Rfam 15.1 is the smaller reference-backed current-family panel.
