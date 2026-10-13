# DAFS DD acceleration combined pilot report（2026-09-19）

## Combined final

`--dd-projected-norm` と `--dd-column-bound` は、どちらも acceleration として不採用、experimental opt-in、default OFF とする。projected-norm は DD iteration `+1.32%`、column-bound は `−0.182%` だったが、column-bound の減少は1 caseの2 mergesだけで、whole/DD wall・CPUおよびCPU19 pinned serialでは全体加速を示さなかった。serialも projected-norm は wall `2.78→2.82 s` / CPU `2.77→2.82 s`、column-bound は wall `2.77→2.87 s` / CPU `2.76→2.86 s` と両方遅くなった。品質・validityは維持されたが、full 680は両戦略とも実行しておらず、full-scale性能やglobal optimalityは主張しない。

freeze後のproduction sourceは変更していない。実行物、OFF control、ON pilot、serial、negative findingを同一reportに保存した。

## 判定

`--dd-projected-norm` は acceleration 候補として不採用とした。pilot の4条件すべてで DD iteration が増え、全体 wall/CPU は改善せず、DD wall は悪化したため、680-run full は実行しない。flag は opt-in の実験用に残し、既定値は変更しない。

candidate の意味は、Polyak denominator から `q_z=0` の boundary outward Z component だけを除外すること。projected update、Q domain、oracle、accuracy filter、stop criteria は変更しない、という main-approved implementation contract で比較した。beam search の出力は近似であり、この結果から global optimality は主張しない。

## 再現性と保存物

| 項目 | 値 |
|---|---|
| old baseline binary | `/tmp/dafs-dd-accel-BpeZy7/dafs-baseline` |
| old baseline SHA-256 | `47ae4d266b60a00291b71aab57abdc2175957d84fa324fc1e2e902170107da1e` |
| historical first-method build path | `build-arm64/src/dafs`（第一pilot時点のpath。現在のsecond binaryを指さない） |
| first-method candidate SHA-256 | `d63ca29621cda1388dbfd3abd62fe7739259562837aa43fcf715053e6e070e67` |
| first-method preserved candidate snapshot | `/tmp/dafs-dd-accel-BpeZy7/dafs-projected-norm` |
| source config | [`luna-relaxed.json`](../benchmarks/configs/luna-relaxed.json) |
| environment | `OMP_NUM_THREADS=1`, seed `42`, ViennaRNA `2.7.2` arm64 |
| beam | all five beam parameters `100` |

既存の all-680 result tree `benchmarks/results-murlet-dataset1-luna-relaxed-arm64-v272` は上書きしていない。今回の result root は `benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272` である。

shared build は `cmake --build build-arm64 -j12` 成功、`ctest --test-dir build-arm64 --output-on-failure` は 10/10 pass。worker 側の focused 4/4、ASAN/UBSAN（`detect_leaks=0`）、USE_ADAGRAD/ADAM rejection、tiny legacy-off regression も pass と報告されている。

## 条件と pilot selection

pilot は `max-iter=1000` の4条件 × 12 cases = 48 paired runs/side。case は12 families と短長のばらつきを固定選択した。

`RF00436-2, RF00031-2, RF00181-2, RF00169-2, RF00059-2, RF00001-2, RF00175-2, RF00050-2, RF00002-2, RF00004-2, RF00230-1, RF00017-1`

| condition | model/sequence | ribosum / final | その他の差分 |
|---|---|---:|---|
| `cl_linear_ath02_iter1000` | `lpc` × LinearAlign | 0.5 / 0.01 | `align-th=0.02` |
| `lpv_probcons_fpct05_iter1000` | `lpv` × LinearProbConsRNA | 0.3 / 0 | `weight=4`, `fold-pct=0.5` |
| `lpc_contralign_ath02_iter1000` | `lpc` × LinearCONTRAlign | 0.3 / 0.005 | `align-pct=0.5`, `align-th=0.02` |
| `lpv_linear_fpct05_iter1000` | `lpv` × LinearAlign | 0.3 / 0 | `align-pct=0.5`, `fold-pct=0.5` |

全条件で `--no-alifold --dynamic-cbp` と align/fold DD beam `100` を固定した。ON config の追加引数は `--dd-projected-norm` だけである。

configs は次の通り。

- [`pilot-baseline`](../benchmarks/configs/luna-dd-acceleration-pilot-baseline.json): old binary, flag OFF
- [`pilot-off-control`](../benchmarks/configs/luna-dd-acceleration-pilot-off-control.json): new binary, flag OFF
- [`pilot-candidate`](../benchmarks/configs/luna-dd-acceleration-pilot-candidate.json): new binary, flag ON
- [`full-baseline`](../benchmarks/configs/luna-dd-acceleration-full-baseline.json) / [`full-candidate`](../benchmarks/configs/luna-dd-acceleration-full-candidate.json): 680-run config（未実行）

## OFF control

old binary と new binary の両方を flag OFF で同じ48 keysに対して比較した。既存 `benchmark_luna_optimization.py compare --mode decoder-only --compare-traces` で time fields を正規化した。

`prediction=48/48`、`quality=48/48`、`summary_non_time=48/48`、`trace_events=48/48`、`errors=0`、`missing=0`。従って、新 binary の既定 OFF path にこの変更由来の legacy regression は観測されなかった。

report: [`comparison-off/luna_optimization_report.json`](../benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/comparison-off/luna_optimization_report.json)

## ON pilot 結果

ON は old baseline → new candidate の順。全48/48が success、metrics valid 48/48。出力変更は failure とせず記録した。prediction と quality dict は 23/48 が byte/dict equal、非時間 summary は 0/48 equal だった。

| condition | DD iter old→ON | MCC Δ | SCI Δ | SPS Δ | whole wall × | DD wall × | CPU × | RSS Δ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `cl_linear_ath02_iter1000` | 7943.6→7961.6 | +0.001687 | +0.013872 | −0.001051 | 1.0128 | 1.0734 | 1.0133 | −0.50% |
| `lpv_probcons_fpct05_iter1000` | 8104.4→8137.2 | −0.000170 | +0.002180 | −0.000321 | 1.0107 | 0.9553 | 1.0111 | −0.44% |
| `lpc_contralign_ath02_iter1000` | 7804.5→8107.8 | +0.000081 | +0.006899 | +0.001582 | 0.9999 | 0.9568 | 0.9999 | −0.69% |
| `lpv_linear_fpct05_iter1000` | 8216.6→8286.8 | −0.001760 | +0.008144 | −0.000846 | 0.9694 | 0.9605 | 0.9704 | +0.30% |

全48 case の単純平均は DD iteration `8017.3→8123.3`（+1.32%）、DD violation `111.65→127.83`、certified violation `275.25→397.50`、certified final gap `10.2915→10.6651`。quality は MCC `0.655921→0.655881`、SCI `0.639855→0.647628`、SPS `0.709746→0.709587` で、大きな quality loss は見えないが、iteration と時間の目的を満たさない。

process resource の全48-run平均は wall `7.7415→7.7675 s`（0.9966×）、DD wall `6.0923→6.1975 s`（0.9830×）、CPU `7.6798→7.7015 s`（0.9972×）、RSS `18360.8→18301.7 KiB`（−0.32%）。`×` は old/candidate であり、1未満は回帰である。

stop reason も `beam_stationary` から `certified_gap`/`max_iterations` 側へ寄り、iteration短縮の証拠にはならなかった。全体結果は [`comparison-on/luna_optimization_report.json`](../benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/comparison-on/luna_optimization_report.json)、条件別CSVは [`luna_optimization_per_condition.csv`](../benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/comparison-on/luna_optimization_per_condition.csv) に保存した。

## CPU19 serial A/B

`RF00169-2 / cl_linear_ath02_iter1000`、metrics-off、CPU19 pinned、baseline/candidate交互順5回で実行した。median は wall `2.78→2.82 s`（0.986×）、CPU `2.77→2.82 s`（0.982×）。全5回が success。ON なので prediction byte diff は全5回で、diff fraction は `0.65787`。短い process の RSS は old `13600 KiB`、candidate median `13408 KiB` だが、memory claim には使わない。

report: [`serial-cpu19/paired.json`](../benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/serial-cpu19/paired.json)

## 第一戦略の再現コマンドと結果扱い

生成は [`benchmark_dd_acceleration.py`](../benchmarks/benchmark_dd_acceleration.py) のみが担当し、実行は既存 runner に委譲した。今回の再現には次を使う。

```sh
python3 benchmarks/run.py \
  --config benchmarks/configs/luna-dd-acceleration-pilot-candidate.json \
  --binary /tmp/dafs-dd-accel-BpeZy7/dafs-projected-norm --jobs 16
python3 benchmarks/benchmark_luna_optimization.py compare \
  --config benchmarks/configs/luna-dd-acceleration-pilot-candidate.json \
  --baseline-results benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/baseline \
  --optimized-results benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/candidate \
  --baseline-kind fresh --mode tight-alignment \
  --report-dir benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/comparison-on
```

`--mode tight-alignment` はこの既存 comparator の「出力変更を non-fatal にする」比較モードとしてだけ使ったもので、config に `--dd-tight-alignment` は付けていない。full 680 は実行していない。第一戦略の結果は棄却済み実験として保存し、defaultは変更しない。

## 第二戦略: `--dd-column-bound`（最終pilot報告、2026-09-19）

`--dd-projected-norm` の棄却は維持する。第二戦略は projected norm を OFF にしたまま、gradient trajectory・row gradient・projected update direction・Q domain を変えず、valid upper bound と certified stop の追加機会だけを調べる `--dd-column-bound` である。raw unordered visitor は sorting せず、full `p_z` + certified-Q coverage、independent column UB、exact `W` の double-outward 評価を使う。raw sparse visitor/column-UB traversal の source-level complexity は `O(L+M)` という限定付きであり、DD iteration、merge、outer refinement、scheduler/IOを含むend-to-end wall-timeの保証ではない。accuracy filter、既存 stop-criteria の閾値、oracle は変更しない。optional refinement は `0` のまま比較した（DD return UB が outer refinement に影響し得るため）。beam search から global optimality は主張しない。

worker の focused 4/4、ASAN/UBSAN bounds+gradient+column smoke、legacy-off tiny exact、CLI invalid-th/nonlinear/both-flags rejection は pass。共有 `cmake --build build-arm64 -j12` は成功し、`ctest --test-dir build-arm64 --output-on-failure` は10/10 passした。第一候補の shared binary `d63ca296...` は第二候補に流用していない。freeze後のproduction sourceは変更していない。

### 再利用する baseline と selection

第一候補で保存済みの old-baseline 48-run tree をそのまま再利用する。baseline の再実行や既存 result の上書きはしない。

- old baseline binary: `/tmp/dafs-dd-accel-BpeZy7/dafs-baseline`
- old baseline SHA-256: `47ae4d266b60a00291b71aab57abdc2175957d84fa324fc1e2e902170107da1e`
- reused baseline results: `benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/baseline`
- frozen column-bound candidate: `/tmp/dafs-dd-column-bound-arm64/dafs-candidate`
- frozen candidate SHA-256: `753129a65753f63b33031c3284ce00be5a09559706c7c21955e6d8d0df2e1695`
- cases: 第一候補と同じ12 cases、4 conditions、計48 paired runs/side
- `max-iter=1000`, beam `100`, `OMP_NUM_THREADS=1`, seed `42`

OFF control は新しい第二候補 binary を flag なしで old baseline と比較する。ON は同じ old baseline に `--dd-column-bound` だけを加える。candidate は output change を許容して記録し、metrics、iterations、stop reasons、violations、gaps、valid metrics、whole/DD wall、CPU、RSS を比較する。

### 使用した config と保存物

生成器 [`benchmark_dd_acceleration.py`](../benchmarks/benchmark_dd_acceleration.py) で、第一候補と同じ12-case selectionの以下のconfigを使用した。pilot のOFF control/ON candidateを実行し、full configはgate FAILのため未実行である。

- [`pilot-baseline`](../benchmarks/configs/luna-dd-column-bound-pilot-baseline.json): 構造確認用。既存48-run baselineを再実行しない
- [`pilot-off-control`](../benchmarks/configs/luna-dd-column-bound-pilot-off-control.json): 新 binary、flag OFF
- [`pilot-candidate`](../benchmarks/configs/luna-dd-column-bound-pilot-candidate.json): 新 binary、flag ON
- [`full-baseline`](../benchmarks/configs/luna-dd-column-bound-full-baseline.json) / [`full-candidate`](../benchmarks/configs/luna-dd-column-bound-full-candidate.json): 680-run config（未実行）
- [`plan`](../benchmarks/configs/luna-dd-column-bound-plan.json): 48-run pilot / 680-run full の固定計画

result root は `benchmarks/results-murlet-dataset1-dd-column-bound-arm64-v272`。preserved candidate snapshot は `/tmp/dafs-dd-column-bound-arm64/dafs-candidate`、SHA-256は `753129a65753f63b33031c3284ce00be5a09559706c7c21955e6d8d0df2e1695` である。old baseline SHAとともに保存した。

### 再現コマンド（実行済み）

new binary の OFF control と ON candidateを既存runnerで実行した。old baseline configは再実行せず、第一候補の結果treeを参照した。

```sh
python3 benchmarks/run.py \
  --config benchmarks/configs/luna-dd-column-bound-pilot-off-control.json --jobs 16
python3 benchmarks/run.py \
  --config benchmarks/configs/luna-dd-column-bound-pilot-candidate.json --jobs 16
python3 benchmarks/benchmark_luna_optimization.py compare \
  --config benchmarks/configs/luna-dd-column-bound-pilot-candidate.json \
  --baseline-results benchmarks/results-murlet-dataset1-dd-acceleration-arm64-v272/pilot/baseline \
  --optimized-results benchmarks/results-murlet-dataset1-dd-column-bound-arm64-v272/pilot/candidate \
  --baseline-kind fresh --mode tight-alignment \
  --report-dir benchmarks/results-murlet-dataset1-dd-column-bound-arm64-v272/pilot/comparison-on
```

OFF control の trace equality は既存 comparator の `--mode decoder-only --compare-traces` で、new binary の flag OFF と reused old baseline を比較する。ON の output difference は failure 扱いせず、quality と validity を別々に判定する。

attribution 用の CPU19 pinned serial は pilot の代表caseで次を使う。runner の baseline は SHA 検証し、candidate-only arg は baseline invocation から除去する。

```sh
python3 benchmarks/benchmark_luna_relaxed.py serial \
  --config benchmarks/configs/luna-dd-column-bound-pilot-candidate.json \
  --baseline-binary /tmp/dafs-dd-accel-BpeZy7/dafs-baseline \
  --candidate-binary /tmp/dafs-dd-column-bound-arm64/dafs-candidate \
  --baseline-sha256 47ae4d266b60a00291b71aab57abdc2175957d84fa324fc1e2e902170107da1e \
  --candidate-only-arg=--dd-column-bound \
  --results-root benchmarks/results-murlet-dataset1-dd-column-bound-arm64-v272/pilot/serial-cpu19 \
  --repetitions 5 --cpu 19
```

### 実行結果（early pilot）

共有build/ctestは上記の通り pass。new binary の flag-OFF control は old baseline と `prediction=48/48`、`quality=48/48`、`summary_non_time=48/48`、`trace_events=48/48`、`errors=0`、`missing=0` で完全一致し、OFF側も success/valid metrics `48/48` だった。ON candidate も全48/48 success、`metrics_validation.valid=48/48`、`quality_complete=48/48`、prediction/quality は48/48 equalだった。非時間 summary は44/48 equalで、差分は上界・gap・証明付き停止に伴う反復数/違反数等だった。

ON の `×` は old/candidate（1より大きいほど candidate が速い）。MCC/SCI/SPS は全条件で差分0だった。

| condition | DD iter old→new | DD violation old→new | certified violation old→new | certified gap old→new | whole wall × | DD wall × | CPU × | RSS Δ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `cl_linear_ath02_iter1000` | 7943.58→7943.58 | 102.083→102.083 | 299.750→299.750 | 13.493525→13.492434 | 1.00480 | 1.00598 | 1.00820 | −0.876% |
| `lpv_probcons_fpct05_iter1000` | 8104.42→8104.42 | 104.917→104.917 | 243.750→243.750 | 8.726607→8.726064 | 1.06352 | 1.02734 | 1.06440 | +0.306% |
| `lpc_contralign_ath02_iter1000` | 7804.50→7746.17 | 106.500→106.583 | 209.333→209.667 | 8.476159→8.435998 | 1.00856 | 0.96959 | 1.00788 | −0.395% |
| `lpv_linear_fpct05_iter1000` | 8216.58→8216.58 | 133.083→133.083 | 348.167→348.167 | 10.469634→10.469634 | 0.91157 | 0.89935 | 0.91155 | −0.276% |

48-run 加重平均は DD iteration `8017.2708→8002.6875`（−14.583/run、−0.182%）。減少は `lpc_contralign_ath02` の `murlet-dataset1-RF00436-2` だけで、merge 5 が `613→82`（−531）、merge 9 が `332→163`（−169）、合計700 iterationを節約した。他11 casesと他3 conditionsは iteration不変。merge 9 の stop reason は `beam_stationary→certified_gap`、merge 5 は `certified_gap` のまま。該当caseの certified gap は merge 5 `0.00704956→0.00466919`、merge 9 `0.02058029→0.00378418` と限定的に縮小した（LPC-CONTRAlign conditionの108 merge summary中、gap数値差分は8件。full pilot全432 merge summaryに対する母数ではない）。

全体resource平均は wall `7.741458→7.813750 s`（0.99075×）、DD wall `6.092325→6.286212 s`（0.96916×）、CPU `7.679792→7.745208 s`（0.99155×）、RSS `18360.8→18304.5 KiB`（−0.307%）。LPV-Linear の DD wall は特に `0.89935×` であり、iteration微減だけから全体加速とは扱えない。

stop reason の変化は LPC-CONTRAlign に限定され、集計は `beam_stationary 20→19`、`certified_gap 1→2`、`max_iterations 87→87`。quality は全48件で byte/dict equal だが、時間目的を満たさない。

CPU19 pinned serial（`RF00169-2 / cl_linear_ath02_iter1000`、metrics-off、交互5回）は全5/5 successかつprediction equal。median は wall `2.77→2.87 s`（0.96516×）、CPU `2.76→2.86 s`（0.96503×）で、controlled timingでも candidate は遅かった。

reports: [`comparison-off`](../benchmarks/results-murlet-dataset1-dd-column-bound-arm64-v272/pilot/comparison-off/luna_optimization_report.json)、[`comparison-on`](../benchmarks/results-murlet-dataset1-dd-column-bound-arm64-v272/pilot/comparison-on/luna_optimization_report.json)、[`serial-cpu19`](../benchmarks/results-murlet-dataset1-dd-column-bound-arm64-v272/pilot/serial-cpu19/paired.json)。

### 最終gate判定

判定は **FAIL**。iteration減少は1 caseの2 mergesに限定され、whole wall/DD wall/CPUは全体で改善せず、CPU19 serialも遅かったため、full 680は実行しない。`--dd-projected-norm` と `--dd-column-bound` はともに experimental opt-in、default OFFのまま。upper-bound theory と certified metadata は保存するが、overall acceleration は実証できなかった。第三候補はこのturnでは追加しない。既存の all-680 tree は変更していない。
