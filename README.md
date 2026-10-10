DAFS: simultaneous aligning and folding of RNA sequences by dual decomposition
==============================================================================

Requirements
------------

* [Vienna RNA package](http://www.tbi.univie.ac.at/~ivo/RNA/) (>= 2.7; also supplies nearest-neighbor energy primitives for linear profile folding)
* (optional)
  [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
  or [Gurobi Optimizer](http://www.gurobi.com/) (>=2.0)
  or [ILOG CPLEX](http://www.ibm.com/software/products/ibmilogcple/) (>=12.0)

Install
-------

	export PKG_CONFIG_PATH=/path/to/viennarna/lib/pkgconfig:$PKG_CONFIG_PATH
	mkdir build && cd build
	cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . 
	cmake --install . # optional

The exact IPknot decoder and `--max-iter=0` require an IP solver. The linear
IPknot path described below works without one.

With `-s lpc` or `-s lpv` together with `--fold-decoder=IPknot` (or
`--ipknot`), DAFS uses a fixed-beam `LinearIPknot` decoder. It keeps IPknot's
threshold layers and selects each layer with `LinearNussinov` over the sparse
LinearPartition support; the support and crossing witnesses are scanned in
linear time for a fixed beam and number of layers. This is the linear folding
decoder used during positive-iteration dual decomposition, so `--max-iter`
only controls the number of coupling iterations. `--max-iter=0` remains the
explicit exact coupled integer-program path and is reported separately.
The linear bound applies to the sparse LinearPartition overload used by DAFS;
the standalone dense matrix overload necessarily has quadratic input size.
With a fixed number of sequences, fixed probability and decoder beams, fixed
positive posterior cutoff, fixed layer count, and fixed DD iteration count,
the sparse decoder and its DAFS multiplier/bound path use linear time and
memory in sequence length. `--dense-lagrangian` and `--max-iter=0` select
separate dense or exact paths. The `--ipknot` option also requests a fixed
number of constrained LinearPartition refolds for the final prediction.

For nonlinear folding models, IPknot continues to use the exact MIP decoder.
The fixed-beam linear decoder is an approximation: it does not provide the
global IP optimum, and the full DAFS run still includes alignment and dual
iteration costs.

For each ungapped input sequence, LinearPartition-C/V caps retained states per
column with its beam and caps internal-loop and helix lengths. Its inside,
outside, and sparse base-pair probability passes therefore use expected
linear time and linear memory in sequence length when those limits are fixed.
The implementation uses hash maps and `nth_element`, so this is not a strict
worst-case timing guarantee. Constrained refolding uses the same sparse
passes; its posterior capacity correction also scans only retained pairs.

Guide tree
----------

The default `--guide-tree=mafft` builds the guide tree inside DAFS for both
linear and nonlinear models. It reimplements MAFFT 7.526's nucleotide
`--globalpair --retree 0 --mixedlinkage 0.1` method: Kimura-200 scoring,
affine global pairwise scores normalized by the smaller self score, and
90% minimum / 10% unweighted-average linkage. MAFFT's terminal-gap rules,
nearest-neighbor tie handling, and left/right child order are preserved.
No MAFFT executable, subprocess, or temporary alignment is required.
Exact pairwise scoring takes quadratic time in sequence length and linear
working memory. For linear-time guide-tree scoring, select `mafft-linear` below.
Both methods store a quadratic distance matrix in sequence count. Guide-tree
time is included in DAFS wall time and the `guide_tree` metrics stage.

* `--guide-tree=mafft` (default): the internally implemented MAFFT-global tree.
* `--guide-tree=mafft-linear --guide-tree-beam=200`: an approximate version of
  the same distance, using affine Viterbi beam search. It keeps at most B
  entries per match/gap state on each antidiagonal, with the same substitution
  scores, opening/extension costs, terminal-gap rules, and mixed linkage as
  exact MAFFT. A fixed B gives expected O(B(L1+L2)) time per sequence pair,
  using O(B+L1+L2) working memory. Only three antidiagonals are retained;
  there is no quadratic fallback. `nth_element` selection gives an expected
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
See `THIRD_PARTY_NOTICES` for MAFFT's BSD license.

With beam 200, all 85 Murlet dataset1 ordered trees and all 3,825 pairwise
distances matched exact scoring. Beam 100 matched 83/85 trees. These are
guide-tree checks; see [validation details](docs/mafft-linear-guide-tree.md)
for the inputs, limitations, and per-input results.

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
