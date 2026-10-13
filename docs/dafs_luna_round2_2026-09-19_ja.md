# DAFS Luna decoder round2 評価報告（2026-09-19）

## 結論

round2 の frozen source から arm64 binary を再構築し、prior-Luna decoder の frozen baseline と decoder-only で比較した。Murlet dataset1 の 680 paired runs（4 linear model combinations × 85 datasets × max-iter 100/1000）で、prediction、quality、非時間 summary、trace event が全件一致した。最終的な max-iter=1000 の wall speedup は cl-linear 1.615×、lpc-ContrAlign 1.589×、lpv-linear 1.457×、lpv-ProbConsRNA 1.535×だった。

最適化 binary は baseline より RSS が約 10.4–12.0% 増加した。これは性能向上と併記すべき実測上のメモリ trade-off であり、特に Beam の長さ 4000 microcase では 34.36 MiB から 59.20 MiB になった。

## 再現条件と identity

- baseline は「original map」ではなく、直前の Luna decoder optimization 後に凍結した prior-Luna decoder の historical result tree: `benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272`。
- baseline binary SHA-256: `49f9abb54daef26701ef0b396192b113c47ee851ddf6c86d23f5d687eee43432`
- optimized binary SHA-256: `75b6fa32805ec02cb3caf32186978cc6f6aed9bf0c51aec410a5347dda9a8531`
- build: `cmake --build build-arm64 -j12`; source-review focused tests and ASAN passed; `ctest --test-dir build-arm64 --output-on-failure`: 9/9 pass（2026-09-19）。
- config: `benchmarks/configs/luna-round2-optimized.json`; `dd_tight_alignment=false`。full run は既存 `benchmarks/run.py` を `--jobs 16` で実行し、各 run は metrics 有効、repetition=1。
- 比較: `benchmark_luna_optimization.py compare --baseline-kind historical --compare-traces --mode decoder-only`。
- new result root: `benchmarks/results-murlet-dataset1-luna-round2-arm64-v272`。既存 baseline tree は上書きしていない。

## Full 680: 出力同一性と性能

comparison の結果は `expected=680, paired=680, errors=0, missing=0, unexpected=0`。prediction/quality/summary_non_time/trace はそれぞれ 680/680 equal だった。MCC、SCI、SPS は baseline と optimized の同値値を示す。時間は baseline → optimized の平均秒、括弧内は baseline/optimized、RSS は KiB と増加率、DD wall は DAFS の DD 部分の実測 wall time である。

| condition | n | MCC / SCI / SPS | wall s (×) | CPU s (×) | RSS KiB (+%) | DD wall s (×) |
|---|---:|---:|---:|---:|---:|---:|
| cl_linear_ath02_iter100 | 85 | 0.685232 / 0.679167 / 0.755282 | 2.9699 → 2.3127 (1.284) | 2.9542 → 2.2971 (1.286) | 16037 → 17890 (+11.55%) | 1.2808 → 0.7452 (1.719) |
| cl_linear_ath02_iter1000 | 85 | 0.691340 / 0.744503 / 0.755934 | 12.5348 → 7.7615 (1.615) | 12.4722 → 7.7014 (1.619) | 16099 → 18036 (+12.03%) | 10.7234 → 6.1357 (1.748) |
| lpv_probcons_fpct05_iter100 | 85 | 0.663473 / 0.660898 / 0.757265 | 3.1306 → 2.3551 (1.329) | 3.1147 → 2.3372 (1.333) | 16583 → 18431 (+11.14%) | 1.3766 → 0.7725 (1.782) |
| lpv_probcons_fpct05_iter1000 | 85 | 0.670270 / 0.725752 / 0.756600 | 13.1225 → 8.5488 (1.535) | 13.0618 → 8.4868 (1.539) | 16625 → 18564 (+11.67%) | 11.2958 → 6.7993 (1.661) |
| lpc_contralign_ath02_iter100 | 85 | 0.685696 / 0.662258 / 0.753633 | 2.9311 → 2.3208 (1.263) | 2.9128 → 2.3048 (1.264) | 16219 → 18058 (+11.34%) | 1.2223 → 0.6904 (1.770) |
| lpc_contralign_ath02_iter1000 | 85 | 0.693255 / 0.718449 / 0.752938 | 11.6051 → 7.3036 (1.589) | 11.5446 → 7.2422 (1.594) | 16268 → 18180 (+11.75%) | 9.7788 → 5.6318 (1.736) |
| lpv_linear_fpct05_iter100 | 85 | 0.669444 / 0.680921 / 0.772789 | 3.3025 → 2.6058 (1.267) | 3.2875 → 2.5888 (1.270) | 17472 → 19291 (+10.41%) | 1.4820 → 0.9501 (1.560) |
| lpv_linear_fpct05_iter1000 | 85 | 0.674727 / 0.764237 / 0.768569 | 14.2131 → 9.7552 (1.457) | 14.1494 → 9.6935 (1.460) | 17574 → 19451 (+10.68%) | 12.4067 → 7.9815 (1.554) |

### Controlled serial A/B

metrics-off、`murlet-dataset1-RF00169-2`、`cl_linear_ath02_iter1000` を baseline/optimized 交互順で 5 回実行した。5 回すべて status=success かつ prediction bytes equal。wall の中央値は 4.63 → 2.89 s（中央値比 1.602×）、CPU の中央値は 4.62 → 2.88 s（1.604×）。RSS は 11840 → 13440 KiB（+13.5%）。

## Length scaling microbenchmark

full run 後、host の heterogeneous core による短時間測定の揺れを避けるため `taskset -c 19` で CPU 19（Cortex-X925、max 4004 MHz）に固定し、warmup=3、repetitions=30、seed=42、beam=100 で再測定した。CPU governor は変更していない。長さは 250/500/1000/2000/4000、入力は O(L) sparse support である。表の speedup は baseline/optimized の wall（CPU もほぼ同値）で、値の順番は L=250,500,1000,2000,4000。

| variant | support mode | raw→unique support（L=250→4000） | wall speedup（L=250→4000） | CPU speedup（L=250→4000） | RSS at L=4000 KiB (old→new) |
|---|---|---|---|---|---:|
| Beam | beam callback | 498→498, 998→998, 1998→1998, 3998→3998, 7998→7998 | 1.642, 1.679, 1.620, 1.662, 1.640 | 1.640, 1.679, 1.620, 1.662, 1.640 | 35184 → 60624 |
| Nuss-certified | preassembled cached support | 681→557, 1369→1120, 2744→2245, 5494→4495, 10994→8995 | 1.306, 1.203, 1.250, 1.344, 1.247 | 1.306, 1.205, 1.250, 1.343, 1.247 | 22200 → 22200 |
| Nuss-uncertified | sparse-matrix wrapper | 679→555, 1367→1118, 2742→2243, 5492→4493, 10992→8993 | 1.167, 1.165, 1.210, 1.251, 1.272 | 1.167, 1.165, 1.216, 1.251, 1.272 | 22224 → 22232 |

15/15 casesで params/input/score/structure/bounds digest がすべて equal、両 binary の output_stable も true だった。certified は preassembled cached support を使い、uncertified は `add_sparse_pair_support` の境界条件を避ける span 3 wrapper とした。star weight は一定値ではなく left に応じて減衰させ、heap/pruning の候補が残る入力にしている。3 repetitions の unpinned 予備測定も `scaling/` に保持するが、短時間測定のため最終値には採用していない。

## 計算量の適用範囲と限界

実装構造から、BeamAlign と LinearNussinov には別々の bound が導ける。BeamAlign は各 layer で最大 3 本の ordered stream から合計 3b 以下の incoming candidate を受け、O(b) merge を行い、b に依存する pruning を行うため、profile length/diagonal count `N=n1+n2` に対して O(N f(b)) となる。LinearNussinov は flat float cache、beam 上の高々 b states に対する suffix `lower_bound`（query あたり binary search O(log(b+1))）、4-pass radix と `C<128` の small-sort、最大 heap size b を使い、候補総数 `C <= O((L+M)b)` から heap ordering を含めて `O((L+M)b log(b+1))` となる。cached-support を事前に組み立てた certified path は `O(L+M)`（固定 b）である。sparse-support core の記法では、固定した正の beam `b > 0`（今回の測定値は b=100）と、重複を含む cached sparse support entry 数 `M` に対して `O(L + M f(b))` の形で表せ、`M=O(L)` の sparse input なら length-linear になる条件である。`std::unordered_map` による score lookup は標準的な expected O(1) とし、callback の計算コストは別に数える。

ここで cached-support とは、support pair を decoder 呼び出し前に重複処理・並べ替え済みの flat/cached 表現へ組み立て、decoder core に渡す条件である。これにより certified 測定は noncached wrapper の前処理を bound に混ぜない。

`Nuss-uncertified` は noncached wrapper の deduplicate/sort を含むため、一般にはその前処理に O(M log M) が加わり得る。したがって wrapper 全体を core bound の証拠とは扱わない。`Nuss-certified` は事前組み立て済み cached-support API を使い、この前処理を切り離して core の測定を行った。

この整理は次を含まない。

- `dense_pair_support` overload は O(L²) の support 構築を行うため対象外。
- BeamAlign の beam=0 unpruned path は対象外。Nussinov 側の beam=0 は `max(1,b)` により effective beam width 1 となる別条件であり、本報告の正の beam bound には含めない。
- fixed positive beam は「b が 100 に固定」という意味ではなく、任意の固定 `b>0` という仮定。今回の実験はその一例として b=100。
- scaling の時間列は実測 evidence であり、漸近 bound の証明ではない。
- whole DAFS や DD iteration count 全体の complexity claim はしていない。full table の DD wall と iteration は観測値としてのみ扱う。

## Artifact index

- [Full new summary](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/summary.csv) and [manifest](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/manifest.json)
- [Full comparison JSON](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/comparison/luna_optimization_report.json) / [CSV](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/comparison/luna_optimization_per_condition.csv)
- [Serial A/B JSON](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/serial-paired/paired.json) / [CSV](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/serial-paired/paired.csv)
- [Pinned scaling JSON](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/scaling-pinned-cpu19-postfull/results/scaling_comparison.json) / [CSV](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/scaling-pinned-cpu19-postfull/results/scaling_comparison.csv)
- [Pinned run identity](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/scaling-pinned-cpu19-postfull/run_context.json)
- [Preliminary unpinned scaling](../benchmarks/results-murlet-dataset1-luna-round2-arm64-v272/scaling/results/scaling_comparison.json)（短時間・異種 core 配置の noise を含むため、最終 timing claim には不使用）
- Reproduction helpers: [round2 harness](../benchmarks/benchmark_luna_round2.py), [scaling source](../benchmarks/luna_round2_scaling.cpp)
