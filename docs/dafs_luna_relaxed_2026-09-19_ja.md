# DAFS Luna relaxed 評価報告（2026-09-19）

## 結論

BeamAlign の選択済み高々 b 個のキーの再ソートと `max_keys_` バッファを削除し、LinearNussinov の `candidate_order` バッファを除去した frozen candidate を、現在の Round2 baseline と比較した。Murlet dataset1 の 8 条件 × 85 cases = 680 paired runs（beam=100、max-iter=100/1000、DD tight=false）で、実行失敗・metrics 不正は 0 件だった。

出力差を failure とは扱わない評価にしたが、今回の実データでは差も 0 件だった。既存の robust exact comparator でも `prediction=680/680`、`quality=680/680`、`summary_non_time=680/680`、`trace_events=680/680` だった。従って quality の MCC/SCI/SPS、DD iteration、violation、bound 関連値は全 paired case で不変である。

性能については、短時間の Beam micro は CPU19 固定で 1.139–1.177× と一貫した。一方、Nussinov certified micro は 0.952–1.034×で、4000 長でも 0.961×だったため、Nussinov 単体の 2–7% 改善とは主張しない。代表的な metrics-off serial AB（CPU19、交互順、5 回）の中央値は wall 2.80→2.77 秒（1.01083×）、CPU 2.80→2.76 秒（1.01449×）であり、これが whole-pipeline に対する主な controlled evidence である。

## identity、build、変更範囲

| 項目 | frozen Round2 baseline | luna-relaxed candidate |
|---|---|---|
| binary SHA-256 | `75b6fa32805ec02cb3caf32186978cc6f6aed9bf0c51aec410a5347dda9a8531` | `47ae4d266b60a00291b71aab57abdc2175957d84fa324fc1e2e902170107da1e` |
| config SHA-256 | `317774efa08fb2aff2976f10a4dfbeefce452322af916dc16579568c0f2d16f1` | `56b16dcabbbe55b03a27bab4d7e3910908268365916610dd96f8b6fe5bb28b53` |
| result tree | `results-murlet-dataset1-luna-round2-arm64-v272` | `results-murlet-dataset1-luna-relaxed-arm64-v272` |

両 tree の selection digest は `faec2c9d1e6644115c6841b4132a0e7fe04c723311315143b1cdae9b36f997dd` で一致し、input SHA-256 と normalized condition args も 680/680 一致した。ViennaRNA は 2.7.2 arm64。`cmake --build build-arm64 -j12` 成功後、`ctest --test-dir build-arm64 --output-on-failure` は 9/9 pass だった。

採用した production change は次の二つだけである。

- BeamAlign: key-sorted layer と stream merge を使い、既存の `nth_element` cutoff/filter を維持したまま、選択済み高々 b 個のキーの再ソートと `max_keys_` バッファを除去した。
- LinearNussinov: candidate の offer-order を保持するだけの `candidate_order` バッファを除去した。探索した nth の追加変更は candidate に含めていない。

baseline/candidate の関連 source manifest は [baseline micro manifest](../benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/baseline-micro-cpu19/source-manifest.sha256) と [candidate micro manifest](../benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/candidate-micro-cpu19/source-manifest.sha256) に固定した。歴史的な Round2 tree は上書きしていない。

## quality、identity、validity

下表の quality は各条件の 85 case 平均。`Δ` は candidate−baseline である。

| condition | MCC | SCI | SPS | Δ(MCC/SCI/SPS) |
|---|---:|---:|---:|---:|
| cl_linear_ath02_iter100 | 0.685232 | 0.679167 | 0.755282 | 0/0/0 |
| cl_linear_ath02_iter1000 | 0.691340 | 0.744503 | 0.755934 | 0/0/0 |
| lpv_probcons_fpct05_iter100 | 0.663473 | 0.660898 | 0.757265 | 0/0/0 |
| lpv_probcons_fpct05_iter1000 | 0.670270 | 0.725752 | 0.756600 | 0/0/0 |
| lpc_contralign_ath02_iter100 | 0.685696 | 0.662258 | 0.753633 | 0/0/0 |
| lpc_contralign_ath02_iter1000 | 0.693255 | 0.718449 | 0.752938 | 0/0/0 |
| lpv_linear_fpct05_iter100 | 0.669444 | 0.680921 | 0.772789 | 0/0/0 |
| lpv_linear_fpct05_iter1000 | 0.674727 | 0.764237 | 0.768569 | 0/0/0 |

- status success: baseline/candidate とも 680/680。
- `metrics_validation.valid`: baseline/candidate とも 680/680。parse/invariant error は 0。
- prediction bytes、quality、非時間 summary、trace events: それぞれ 680/680 equal。
- trace 平均は `dd_iterations` 4570.93→4570.93、`dd_violations_total` 202.87→202.87、`certified_violations_total` 304.91→304.91。後者は全 iteration の累積ではなく、従来定義どおり merge ごとの最後の iteration を合計した値である。
- 追加 micro validator は 15/15 で stable、score/structure/bounds digest equal。Nussinov の finite bound と `score ≤ upper ≤ {pruned, additive}` の検査は 10/10 pass。full run は runner の metrics validity と exact trace 比較で検証した。

## full 680 の timing と memory

full run は jobs=16 の historical baseline tree 対 candidate tree であり、host の異種 CPU・実行時 contention を含む observational timing である。下表の倍率は「平均 baseline 秒 / 平均 candidate 秒」で、case ごとの倍率の平均ではない。

| condition | mean wall 秒（baseline→candidate、×） | mean CPU 秒（baseline→candidate、×） | max RSS KiB（baseline→candidate、増加率） |
|---|---:|---:|---:|
| cl_linear_ath02_iter100 | 2.3127→2.3376 (0.989) | 2.2971→2.3205 (0.990) | 17890→18004 (+0.64%) |
| cl_linear_ath02_iter1000 | 7.7615→7.4360 (1.044) | 7.7014→7.3772 (1.044) | 18036→18162 (+0.70%) |
| lpv_probcons_fpct05_iter100 | 2.3551→2.4153 (0.975) | 2.3372→2.3992 (0.974) | 18431→18544 (+0.62%) |
| lpv_probcons_fpct05_iter1000 | 8.5488→8.2907 (1.031) | 8.4868→8.2265 (1.032) | 18564→18668 (+0.56%) |
| lpc_contralign_ath02_iter100 | 2.3208→2.4072 (0.964) | 2.3048→2.3882 (0.965) | 18058→18154 (+0.53%) |
| lpc_contralign_ath02_iter1000 | 7.3036→7.0127 (1.041) | 7.2422→6.9519 (1.042) | 18180→18249 (+0.38%) |
| lpv_linear_fpct05_iter100 | 2.6058→2.5024 (1.041) | 2.5888→2.4860 (1.041) | 19291→19438 (+0.76%) |
| lpv_linear_fpct05_iter1000 | 9.7552→9.2231 (1.058) | 9.6935→9.1607 (1.058) | 19451→19573 (+0.63%) |

full RSS は約 +0.4–0.8% の観測値であり、保存メモリ量の保証ではない。CPU19 micro は baseline 24296 KiB、candidate 25532 KiB という別プロセスの parent RSS floor の影響を受け得る観測であり、その差 1236 KiB を実装固有の増加に帰属できない。

## controlled micro と serial AB

micro は `taskset -c 19`、Cortex-X925、warmup=3、repetitions=30、seed=42、beam=100。値の順は length=250/500/1000/2000/4000 の wall speedup である。

| variant | wall speedup | 判定 |
|---|---|---|
| Beam | 1.139 / 1.177 / 1.157 / 1.158 / 1.167 | robust な cutoff/filter gain |
| Nuss-certified | 0.952 / 1.021 / 1.034 / 1.034 / 0.961 | mixed、overall gain は主張しない |
| Nuss-uncertified | 0.997 / 0.986 / 1.014 / 1.014 / 1.029 | mixed、補助観測 |

serial AB は metrics-off、`cl_linear_ath02_iter1000`、`murlet-dataset1-RF00169-2`、CPU19、baseline/candidate の交互順で 5 回。全 5 回が success、prediction byte diff は 0。wall speedup は 1.007–1.018×、CPU speedup は 1.007–1.014×、中央値は上記の 1.01083× / 1.01449×。candidate RSS の中央値は 13600 KiB、baseline は 13248 KiBだが、短時間 process の RSS floor を含むため memory claim には用いない。

## 計算量の条件と限界

Beam の各 layer は高々 3b 個の候補を扱い、選択済み高々 b 個のキーの再ソートを `max_keys_` バッファなしで処理する。選択処理のコストを f(b) と置けば、`nth_element` の worst-case を含めても selection cost は length に依存せず、N=n1+n2 に対して従来の bound `O(N f(b))`、memory `O(Nb)` を維持する。ここで f(b) は b に依存する選択・filter・merge のコストであり、nth の平均ケースや「十分な候補数」を前提にしない。Nussinov の変更は offer-order の保持を省く定数因子削減で、従来の `O((L+M)b log(b+1))`（M は sparse support entry 数、hash lookup は expected O(1)）を変えない。固定 b かつ M=O(L) の sparse support では length-linear になる。noncached support の dedup/sort には追加コストがあり得るため core bound から分離し、dense support は対象外とする。Beam の beam=0 は unpruned path、Nussinov の beam=0 は effective beam=1 という別条件であり、whole DAFS の漸近計算量や DD iteration の短縮は主張しない。

## 最小の評価適応と再現手順

小さな exact/unit test は feasibility、traceback、certificate、bound invariant の検査として維持する。end-to-end の候補評価では prediction identity を必須 failure にせず、次を必須 gate とする。

- baseline/candidate の result set が非空で、期待件数 680、selection/config/input/condition identity を検証する。
- runner validator で status、parse/invariant error、finite bound と validity を fail-closed に検証する。
- MCC/SCI/SPS、prediction byte diff count/fraction、DD iteration、violation、gap を paired で報告する。
- exact identity が必要な回帰確認は既存 `benchmark_luna_optimization.py --mode decoder-only --compare-traces` を補助的に実行する。

再現コマンド:

```bash
cmake --build build-arm64 -j12
ctest --test-dir build-arm64 --output-on-failure
python3 benchmarks/run.py --config benchmarks/configs/luna-relaxed.json --jobs 16
python3 benchmarks/benchmark_luna_relaxed.py compare \
  --baseline-results benchmarks/results-murlet-dataset1-luna-round2-arm64-v272 \
  --candidate-results benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272 \
  --output-root benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/comparison-approx \
  --expected-runs 680
```

主要 artifacts: [luna-relaxed config](../benchmarks/configs/luna-relaxed.json)、[full summary](../benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/summary.csv)、[candidate manifest](../benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/manifest.json)、[approx comparison JSON](../benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/comparison-approx/comparison.json)、[exact supplemental report](../benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/comparison-exact-supplemental/luna_optimization_report.json)、[serial AB](../benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/serial-paired-cpu19/paired.json)、[CPU19 micro comparison](../benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272/micro-comparison-cpu19/micro_comparison.csv)、[new helper](../benchmarks/benchmark_luna_relaxed.py)。
