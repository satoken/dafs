DAFS: simultaneous aligning and folding of RNA sequences by dual decomposition
==============================================================================

Requirements
------------

* [Vienna RNA package](http://www.tbi.univie.ac.at/~ivo/RNA/) (>= 2.7; also supplies nearest-neighbor energy primitives for linear profile folding)
* (optional)
  [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
  or [Gurobi Optimizer](http://www.gurobi.com/) (>=2.0)
  or [ILOG CPLEX](http://www.ibm.com/software/products/ibmilogcple/) (>=12.0),
  SCIP, or HiGHS

Install
-------

	export PKG_CONFIG_PATH=/path/to/viennarna/lib/pkgconfig:$PKG_CONFIG_PATH
	mkdir build && cd build
	cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . 
	cmake --install . # optional

CMake builds the bundled IPknot decoder in `vendor/ipknot` by default; no
separate IPknot checkout, network download, or CMake path is required. The
snapshot and per-file hashes are recorded in `vendor/ipknot/UPSTREAM.json`.
To override the bundled snapshot with another checkout:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DIPKNOT_SOURCE_DIR=/path/to/ipknot
cmake --build build -j
ctest --test-dir build --output-on-failure
```

For an installed decoder, set `-DIPKNOT_SOURCE_DIR=` and
`-DCMAKE_PREFIX_PATH=/path/to/ipknot-decoder-install` instead. The decoder's
public headers and exported CMake target are sufficient for this mode.

Only the explicit exact coupled `--max-iter=0` path requires an IP solver.
Nussinov and the layered IPknot DD decoder work without an IP solver for all
probability models, including CONTRAfold. Solver
selection uses `ENABLE_CPLEX`, `ENABLE_GUROBI`, `ENABLE_SCIP`, or `ENABLE_HIGHS`;
otherwise GLPK is discovered automatically. These settings apply to both DAFS
and the decoder built from source. An installed decoder retains its own
build-time solver selection. Set `-DENABLE_ILP=OFF` for a solver-free source
build; the IPknot DD path described below works without a solver.

With `--fold-decoder=IPknot` (or `--ipknot`), DAFS uses the reusable
IPknot decoder's bounded beam dual
decomposition backend regardless of the probability model (`-s CONTRAfold`,
`-s lpc`, `-s lpv`, etc.), including auxiliary BPP inputs. This changes decoding
only; CONTRAfold still calculates the same base-pair probabilities.
`--fold-dd-beam` and `--fold-final-beam` set its
Nussinov beam for coupling and final prediction, respectively. Its other DD
settings use upstream defaults: 50 inner iterations, crossing beam 100, and
16 witnesses. These inner iterations are independent of DAFS's `--max-iter`,
which controls the outer alignment/folding coupling iterations.
`--max-iter=0` remains the explicit exact coupled integer-program path and is
reported separately. The `--ipknot` option also requests a fixed number of
constrained refolds with the selected probability model for the final prediction.

DD is an approximate pseudoknot decoder; its prediction can differ from the
previous nonlinear ILP backend. DAFS uses relaxed upper bounds for the outer
coupling problem and reports uncertified agreement as `beam_stationary`.
For example, no additional decoder option or solver is required:

```sh
dafs -a CONTRAlign -s CONTRAfold --ipknot --refinement=0 input.fa
```

The adapter receives DAFS's objective `w * (posterior - threshold) - multiplier`,
including pairs present only in the multiplier support. DAFS preserves summed
level scores rather than IPknot's default level average. Returned symmetric
partners and level assignments are converted to DAFS's one-sided structures
and bracket strings.

All structure decoders allow isolated base pairs by default during progressive,
final, and constrained refinement decoding, including the exact coupled
`--max-iter=0` path. Use `--disallow-isolated` to require every selected pair to
have an immediately adjacent stacked pair. The older `--allow-isolated` option
remains accepted for compatibility, and `--allow-isolated=false` also disables
isolated pairs. If both options are present, `--disallow-isolated` takes
precedence. The constrained recurrences retain below-threshold candidates when
they support a profitable helix. Missing posterior and multiplier entries are
not added as artificial stack partners. LinearNussinov keeps paired states on
the sparse candidate support alongside its interval beam, preserving linear
work and storage for fixed beams and sparse degrees. DD repair removes pairs
that lose stack support on either aligned profile.

For nonlinear folding models, the default `DAFSNussinov` uses DAFS's original
candidate-list dynamic program, restored from `d7f4b07` before IPknot kernel
integration. Its interval scores, addition order and tie precedence retain the
previous float-DP behavior. Only positive local scores enter the candidate
lists when isolated pairs are allowed. Sparse posterior inputs use the
historical dense interval chart, and sparse profile bonuses are added before
the threshold is subtracted. The stacked-pair recurrence continues to enforce
`--disallow-isolated`, including below-threshold helix partners. Both modes
use the original serial DAFS kernel, with cubic worst-case time and quadratic
chart memory. Run-start metrics identify this backend as `nussinov_backend=DAFS`.

Finite-beam `LinearNussinov` always uses the bundled IPknot `DDNussinov`
improved-beam kernel, including short and very sparse inputs. Suffix dominance
runs before finite-beam selection; exact closed-component sparsification is
restricted to unpruned decoding. DAFS's minimum pair span, candidate support,
tie order and isolated-pair policy are preserved by the adapter. When all inside
states have been pruned, native improved beam may still form a pair with an
empty interior; finite-beam structures and the outer DD trajectory can differ.
IPknot's pruning observer feeds DAFS's existing completion bounds, so the
outer DD termination checks remain certified. Pair weights and masks are
refreshed on each solve; `reset_pairs()` reuses the chart and traceback buffers
when endpoints change at the same length. Run-start metrics identify the
linear backend as
`nussinov_backend=IPknot`. The independent dense `Nussinov` and the IPknot
exact sparse adapter `SparseNussinov` remain direct test/reference APIs.

These defaults use `--fold-decoder=Nussinov` (also the CLI default), with no
extra runtime option or ILP solver needed. `--fold-decoder=IPknot` selects the
separate layered pseudoknot decoder described above.
These choices apply to decoding independently of BPP calculation.

`--fold-dp-threads N` selects within-DP workers for progressive, recovery and
final LinearNussinov decoding, and for IPknot DD. The restored nonlinear
DAFSNussinov kernel stays serial; `--dd-threads` can distribute its independent
DD subproblems instead. The default thread request is 1; 0 uses the
OpenMP runtime limit (for example, `OMP_NUM_THREADS=4`). Negative values are
rejected. Builds without OpenMP accept only 0 and 1. The option does not change
ILP solver threads or probability calculation. The direct IPknot exact
reference API parallelizes interval charts by interval length; very sparse
exact charts remain serial to avoid allocating quadratic interval storage.
LinearNussinov keeps one team across endpoints,
uses private candidate buffers and merges ties in original support order.
Like IPknot, it falls back to serial execution below 8,192 estimated prefix
offers per endpoint and caps workers at one per 4,096 offers. Small/sparse
inputs therefore need not use all requested workers or become faster.
For fixed beams and workers, sparse beam storage remains O(L*b + M + N*L).
Scores, structures and DAFS beam certificates are deterministic across worker
counts; the independent dense Nussinov API also supports interval parallelism.
No global OpenMP settings or nested worker teams are introduced.

`--dd-threads N` instead distributes each DD iteration's X folding, Y folding
and Z alignment subproblems across up to three workers. The default is 1; 0
uses the OpenMP runtime limit. Each folding subproblem owns its decoder and
workspace. Objectives are summed in X/Y/Z order after joining, before shared
multiplier updates, certificates or primal recovery. Results remain
deterministic across worker counts. If `--fold-dp-threads` is anything other
than 1, within-DP parallelism takes precedence and DD subproblems run serially.
Without OpenMP, only 0 and 1 are accepted. Metrics record the requested
`dd_threads` at run start and the actual maximum `dd_workers` for each DD merge.
For example, `dafs --dd-threads 3 input.fa` assigns one worker per subproblem;
it does not parallelize probability calculation, final folding or profile merges.
Pair-match variables are evaluated after dynamic CBP generation and pricing,
which depend on the X/Y/Z solutions.

An external IPknot checkout/package must expose `DDNussinov::PruningObserver`,
and `reset_pairs()`. The bundled snapshot already
includes these interfaces.

For IPknot, the isolated-pair setting applies to both DD and ILP. With
`--disallow-isolated`, DD requires an immediately adjacent stacked pair on the
same level; ILP requires neighboring paired left endpoints and neighboring
paired right endpoints on that level. The DD
backend can produce different predictions from DAFS's previous greedy
layer-by-layer approximation. It is approximate and does not guarantee the
full ILP optimum. Its bounds refer to its retained crossing/witness graph.

For fixed levels, beams, witness widths, inner iteration counts, and recovery
budgets, bounded DD uses expected linear time and memory in sequence length
plus retained sparse candidates. The DAFS adapter scans sparse probability
and multiplier supports without densifying them. The dense overload has
quadratic input size. With a fixed number of sequences, fixed probability
beams, positive posterior cutoff, and fixed outer DD iterations, the sparse
decoder and DAFS multiplier/bound path remain linear in sequence length.
`--dense-lagrangian` and `--max-iter=0` select separate dense or exact paths;
the full DAFS run also includes guide-tree and alignment costs.

For each ungapped input sequence, LinearPartition-C/V caps retained states per
column with its beam and caps internal-loop and helix lengths. Its inside,
outside, and sparse base-pair probability passes therefore use expected
linear time and linear memory in sequence length when those limits are fixed.
The implementation uses hash maps and `nth_element`, so this is not a strict
worst-case timing guarantee. Constrained refolding uses the same sparse
passes; its posterior capacity correction also scans only retained pairs.

`-a LinearAlign-CONTRAlign` retains CONTRAlign's complete five-state,
double-affine RNA model and all 24 learned parameters. Its forward/backward
passes use the same retained cell graph, with parameters loaded from the
bundled `contralign/Defaults.ipp`. Positive `--align-beam` controls pruning;
the `linearcontralign5_parity` test checks unpruned posteriors against the
nonlinear CONTRAlign implementation, including both affine gap families.
`--metrics-jsonl` records `alignment_probability_states` (5 for CONTRAlign).

Guide tree
----------

The default `--guide-tree=mafft-dafs` keeps DAFS's posterior-derived sequence
distances, `d(i,j) = 1 - sim(i,j)`, and applies MAFFT's mixed linkage: 90%
minimum plus 10% unweighted average. It uses the same linkage routine, nearest
neighbor tie handling, and left/right child order as the MAFFT-global mode.
The similarity matrix used by PCT and scoring is unchanged. This works for
both linear and nonlinear models without an additional global pairwise pass.
`--guide-tree=dafs` selects the original DAFS linkage.

The explicit `--guide-tree=mafft` mode reimplements MAFFT 7.526's nucleotide
`--globalpair --retree 0 --mixedlinkage 0.1` method: Kimura-200 scoring,
affine global pairwise scores normalized by the smaller self score, and
90% minimum / 10% unweighted-average linkage. MAFFT's terminal-gap rules,
nearest-neighbor tie handling, and left/right child order are preserved.
No MAFFT executable, subprocess, or temporary alignment is required.
Exact pairwise scoring takes quadratic time in sequence length and linear
working memory. For linear-time guide-tree scoring, select `mafft-linear` below.
Both methods store a quadratic distance matrix in sequence count. Guide-tree
time is included in DAFS wall time and the `guide_tree` metrics stage.

* `--guide-tree=mafft-dafs` (default): DAFS distances with MAFFT mixed linkage.
* `--guide-tree=mafft`: the internally implemented MAFFT-global tree.
* `--guide-tree=mafft-linear --guide-tree-beam=200`: an approximate version of
  the same distance, using affine Viterbi beam search. It keeps at most B
  entries per match/gap state on each antidiagonal, with the same substitution
  scores, opening/extension costs, terminal-gap rules, and mixed linkage as
  exact MAFFT. A fixed B gives expected O(B(L1+L2)) time per sequence pair,
  using O(B+L1+L2) working memory. It first tries the row-wise MAFFT kernel
  in a band with half-width at most 16B, using at most two passes. A bound
  on every excluded path must certify the exact score before that result
  is used. Otherwise it runs the beam search. The band budget is fixed in B;
  it does not expand without limit on long inputs. Beam fallback retains
  three antidiagonals, uses coordinate-indexed scores, and selects entries
  in place. `nth_element` selection gives an expected
  rather than strict worst-case time bound. Ties obey the hard beam cap.
  Larger positive beams retain more candidates and usually approach exact
  scoring (not necessarily monotonically); B above the shorter sequence
  length disables pruning and recovers the exact score. Pair orientation is
  canonicalized so reversing the two sequences gives the same distance.
  Beam size is independent of `--align-beam` and `--linfold-beam` used by the
  probability models.
* `--guide-tree=dafs`: the original posterior-similarity tree.
* `--guide-tree=raf`: the internally implemented RAF 1.00 tree. Distances use
  raw alignment posteriors before PCT, normalized by the shorter sequence.
  This compatibility mode intentionally retains RAF 1.00's upper-triangle-only
  initialization and subsequent lower-triangle reads, which can produce
  input-order-dependent chain trees. It matches the RAF mode benchmarked here.

All methods align the original input sequences and preserve their names and
input order. The guide-tree choice does not replace the similarities used to
weight PCT. `--metrics-jsonl` records the selected `guide_tree_method` and
`guide_tree_beam` (0 for the exact/non-beam methods).
`--guide-tree-diagnostics` prints the full-precision DAFS distance and
similarity for every input pair to stderr, for matched-distance comparisons.
`--raw-align-cutoff` optionally controls raw CONTRAlign/LinearAlign posterior
pruning independently of the decoder's `--align-th`. Its default, -1, follows
`--align-th`; profile matching probabilities retain their fixed 0.01 cutoff.
The benchmark recommendations use `--raw-align-cutoff=0.01 --align-th=0.001`.
See `THIRD_PARTY_NOTICES` for MAFFT's BSD license.

With beam 200, all 85 Murlet dataset1 ordered trees and all 3,825 pairwise
distances matched exact scoring. The optimized implementation also matched
the 16S/23S examples, including all 20 pairwise distances. On those two inputs
it was 8.63x/3.66x faster than the initial beam implementation; these are
guide-tree measurements, not complete DAFS runtimes. See
[validation details](docs/mafft-linear-guide-tree.md) for the inputs,
limitations, and per-input results.

Usage
-----
	DAFS: dual decomposition for simultaneous aligning and folding RNA sequences.
	Usage:
  	  dafs [OPTION...] FILE
    	
  	  -h, --help             Print usage
	  -w, --weight arg       Weight of the expected accuracy score for
	                         secondary structures (default: 4.0)
	      --ribosum-weight arg
	                         Weight of RIBOSUM85-60 pair-pair match scores
	                         (default: 0.075; use 0 to disable)
	      --final-ribosum-weight arg
	                         Weight of the RIBOSUM85-60 self-profile bonus in
	                         final Nussinov decoding (default: 0)
  	  -m, --max-iter T       The maximum number of iteration of the subgradient 
                         	  optimization (default: 600)
  	  -v, --verbose arg      The level of verbose outputs (default: 0)
    
 	  Aligning options:
	  -a, --align-model arg     Alignment model for calculating matching
	                           probabilities (value=CONTRAlign, ProbCons,
	                           LinearAlign, LinearAlign-CONTRAlign,
	                           LinearAlign-ProbCons,
	                           LinearAlign-ProbConsRNA)
	                           (default: ProbCons)
  	  -u, --align-th arg        Threshold for matching probabilities (default: 
                            	  0.01)
    
 	  Folding options:
  	  -s, --fold-model arg     Folding model for calculating base-pairing 
                              probablities (value=Boltzmann, Vienna,
                              CONTRAfold, lpv, lpc) (default: Boltzmann)
  	  -t, --fold-th arg        Threshold for base-pairing probabilities 
                           	  (default: 0.2)
      -T, --fold-th1 arg       Threshold for base-pairing probabilities of the
                              conclusive common secondary structures
          --alifold           Use profile base-pairing probabilities
                              (LinearAlifold for lpv/lpc, RNAalifold otherwise;
                              disabled by default)
          --no-alifold        Disable profile base-pairing probabilities
                              (default)
          --alifold-stages arg
                              Select profile folding stages: none, progressive,
                              final, or both. This is mutually exclusive with
                              --alifold/--no-alifold; default behavior is both
                              when --alifold is given.
          --allow-isolated     Allow isolated base pairs (default; compatibility option)
          --disallow-isolated  Prohibit isolated base pairs in structure decoding
      	  --ipknot             Set optimized parameters for IPknot decoding 
                           	  (--fold-decoder=IPknot -g4,8 -G2,4 --bp-update1)

`--final-ribosum-weight rho` changes the final Nussinov pair score to
`p(i,j) - threshold + rho * R(i,j)`, where `R(i,j)` is the expected
RIBOSUM85-60 score between two draws from the aligned pair-type distribution.
Gaps and ambiguous residues contribute zero while the denominator remains the
total number of aligned sequences, so gappy column pairs are attenuated
quadratically.  The bonus is evaluated on the sparse BPP support and is disabled
by default (`rho=0`).

Experimental linear-DD controls (all disabled by default):

- `--dd-recovery-interval=10` proposes alignments from an eight-iterate
  multiplier window and from the certified track every ten iterations, then
  verifies and scores repaired feasible solutions using the original objective.
- `--dd-beam-projected-norm` changes the beam track's norm only;
  `--dd-beam-eta=0.75` overrides its Polyak scale without changing the certified
  track's scale. These require linear alignment and folding DD decoders.
- `--dd-diagnostics --metrics-jsonl=FILE` serializes fixed merges of length at
  most nine for offline exhaustive/LP diagnosis; no exact oracle is added to
  the production linear algorithm.
- `--dd-unpruned-bound` uses the full max-alignment value as a certificate
  only when that particular beam call discarded no states.
- `--dd-block-bound=32` adds fixed-width row-block monotone alignment upper
  bounds on both multiplier tracks (width 1–64; 0 disables it). Inter-block
  constraints are relaxed; fixed-width work is linear in sparse support size.
- `--dd-outward-lb` verifies and downward-rescores the best feasible solution.
  This addresses lower-bound rounding stalls, but can change subsequent
  Polyak updates and predictions. It does not interval-certify every existing
  upper-bound/CBP arithmetic path.

These options do not loosen certificate tolerances or guarantee faster
convergence. They are experimental controls and may alter runtime or the
resulting alignment; compare quality against the fixed baseline before using
them in production.

Example
-------

	% dafs RF00005:0.fa
	[ 0.0985233 [ 0.585795 [ 0.933469 M68929-1/151018-150946 X00360-1/1-73 ] [ 0.826623 X12857-1/421-494 [ 0.935672 J05395-1/2325-2252 M16863-1/21-94 ] ] ] [ 0.349897 [ 0.780743 J04815-1/3159-3231 [ 0.96716 J01390-1/6861-6932 M20972-1/1-72 ] ] [ 0.74278 K00228-1/1-82 AC009395-7/99012-98941 ] ] ]
	>SS_cons
	(((((((...(((..............))).......(((((..........)))))......(.((((.......))))).))))))).
	> J01390-1/6861-6932
	CAGGUUA-GAGCC-AGGU-GGU-UA--GGCG------UCUUGUUU--GG-GUCAAGA-AAUUGU-UAUGUUCGAAUCAUAA-UAACCUGA
	> J05395-1/2325-2252
	GGUUUCG-UGGUC-UAGUCGGUUAU--GGCA------UCUGCUUA--AC-ACGCAGA-ACGUCC-CCAGUUCGAUCCUGGG-CGAAAUCG
	> K00228-1/1-82
	GGUUGUUUG-GCCGA-GC-GGU-CUAAGGCGCCUGAUUCAAGCUCAGGU-AUCGUAA--GAUGCAAGAGUUCGAAUCUCUU-AGCAACCA
	> AC009395-7/99012-98941
	GGCUCAA-U----------GGU-CUAG-GGGUAUGAUUCUCGCUUUGGG-UGCGAGA--GGUCC-CGGGUUCAAAUCCCGG-UUGAGCCC
	> J04815-1/3159-3231
	AGAGCUU-GCUCC-CAAA-GCU-UG--GGUG------UCUAGCUG--AU-AAUUAGA-CUAUCA-AGGGUUAAAUUCCCUUCAAGCUCUA
	> M20972-1/1-72
	AGGGCUA-UAGCU-CAGC-GGU-AG--AGCG------CCUCGUUU--AC-ACCGAGA-AUGUCU-ACGGUUCAAAUCCGUA-UAGCCCUA
	> M68929-1/151018-150946
	CGCGGGA-UAGAG-UAAUUGGU-AA--CUCG------UCAGGCUC--AU-AAUCUGA-AUGUUG-UGGGUUCGAAUCCGAC-UCCCGCCA
	> X00360-1/1-73
	CCGACCU-UAGCU-CAGUUGGU-AG--AGCG------GAGGACUG---UAGAUCCUU-AGGUCA-CUGGUUCGAAUCCGGU-AGGUCGGA
	> X12857-1/421-494
	GCGGAUG-UAGCC-AAGUGGAUCAA--GGCA------GUGGAUUG--UG-AAUCCACCAUG-CG-CGGGUUCAAUUCCCGU-CAUUCGCC
	> M16863-1/21-94
	GGGCUCG-UAGCU-CAGAGGAUUAG--AGCA------CGCGGCUA--CG-AACCACG-GUGUCG-GGGGUUCGAAUCCCUC-CUCGCCCA


References
----------

* Sato, K., Kato, Y., Akutsu, T., Asai, K., Sakakibara, Y.: DAFS: simultaneous aligning and folding RNA sequences via dual decomposition. *Bioinformatics*, 28(24):3218-3224, 2012.
* Klein, R. J., Eddy, S. R.: RSEARCH: finding homologs of single structured RNA sequences. *BMC Bioinformatics*, 4:44, 2003.
