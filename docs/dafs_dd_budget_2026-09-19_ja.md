# DD maxiter 予算・時間対品質評価（2026-09-19）

## 結論

現行バイナリで `maxiter=100` が `baseline=500` と同等になる、または全条件で
`5x` の品質対時間改善になる、という証拠は得られなかった。全85ケースの
baseline 100 対 baseline 500 は、平均で `ΔMCC=-0.005236`、`ΔSCI=-0.056214`、
`ΔSPS=+0.001785` であり、SCI が明確な識別指標になった。

パイロット後の条件別フォローアップでは、次の2つを記録対象とした。これは
全体ゲートを通った普遍的な選択ではなく、12ケースを含むパイロットで見えた
条件別の候補である。

* `projected-norm, iter=200, lpc_contralign_ath02`
* `projected-norm, iter=300, lpv_probcons_fpct05`

`projected-norm` は現時点でも opt-in の候補に留め、デフォルトは変更しない。
実運用での採用判断は、後述の metrics-off 固定CPU測定と、内部DD収束指標の
許容範囲を合わせて行うべきである。以下の閾値は絶対平均差 `0.005`（参考として
`0.001`）を表示するための記述的な基準であり、同値性検定、統計的有意性、
生物学的同等性、または未観測データへの一般化を主張しない。

## 実験範囲と再現条件

* データ: Murlet dataset1 の85ケース（17 family × 5）。パイロットは固定12ケース。
* モード: baseline、`--dd-projected-norm`、`--dd-column-bound`、
  `--dd-tight-alignment` を個別評価。
* 予算: `maxiter=100, 200, 300, 500`。maxiter は merge ごとの上限であり、
  100 対500は概ね DD反復上限900対4500相当として解釈する（実際の
  `dd_iterations` はケースごとに記録）。
* 固定値: beams=100、seed=42、OMP_NUM_THREADS=1、refinement=0、Vienna=2.7.2。
* バイナリ: `build-arm64/src/dafs`、SHA256
  `753129a65753f63b33031c3284ce00be5a09559706c7c21955e6d8d0df2e1695`。
* 実行数: パイロット768、全85 baseline 100/500 が680、条件別フォローアップが
  340、合計1,788 run。全5440 runの無条件実行はしていない。

比較は同一 dataset・condition・repetition の paired delta（candidate−reference）
である。平均、中央値、p10/p90、最小/最大、符号数、個別 worst case を
`dd_budget_report.json` に保存した。

## 主要2条件の全85平均

表の `mean` は各 arm の実測平均、`Δ` は左の projected-norm arm から右の
reference arm を引いた値である。`viol` は内部DD consensus violation の
`dd_violations_total`、`gap` は最後の progressive merge の certified `bestUB−LB`
（`dd_final_gap_last`）、
`DDwall` は `dd_wall_seconds`、`wall` は全 wall time である。これらは品質指標と
別の収束・コスト指標である。

### LPC×LinearCONTRAlign / projected-norm 200

| arm | MCC | SCI | SPS | viol | gap | DDwall (s) | wall (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| baseline 500 | 0.693145 | 0.711410 | 0.751998 | 108.059 | 5.706 | 2.771 | 4.399 |
| baseline 200 | 0.689518 | 0.692110 | 0.752334 | 167.882 | 7.517 | 1.101 | 2.520 |
| projected-norm 200 | 0.690427 | 0.706440 | 0.751533 | 162.776 | 7.950 | 1.096 | 2.516 |
| projected 200 − baseline 500 | −0.002718 | −0.004970 | −0.000465 | +54.718 | +2.244 | −1.675 | −1.883 |
| projected 200 − baseline 200 | +0.000909 | +0.014330 | −0.000802 | −5.106 | +0.433 | −0.005 | −0.004 |

対 baseline500 の SCI は緩い記述基準内だが、同時に violation と gap は大きい。
対 same-budget baseline では SCI は改善する一方、gap は改善していない。

### LPV×LinearProbCons / projected-norm 300

| arm | MCC | SCI | SPS | viol | gap | DDwall (s) | wall (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| baseline 500 | 0.670141 | 0.717345 | 0.756692 | 127.412 | 5.005 | 3.351 | 4.880 |
| baseline 300 | 0.668758 | 0.704590 | 0.756037 | 173.471 | 5.991 | 1.993 | 3.488 |
| projected-norm 300 | 0.668864 | 0.715636 | 0.756061 | 163.141 | 6.500 | 2.076 | 3.553 |
| projected 300 − baseline 500 | −0.001277 | −0.001708 | −0.000631 | +35.729 | +1.494 | −1.275 | −1.326 |
| projected 300 − baseline 300 | +0.000106 | +0.011046 | +0.000025 | −10.329 | +0.508 | +0.084 | +0.065 |

LPV の projected-norm 300 は nonpilot73 の SCI 差も −0.001661 と緩い基準内だが、
全85平均の gap は baseline500 より大きい。品質の近さを DD consensus の近さと
取り違えないことが重要である。

## パイロット・層別結果

パイロットで得た selected condition の baseline500 対差は、LPC が
`ΔMCC=-0.004018, ΔSCI=-0.000061, ΔSPS=+0.001612`、LPV が
`ΔMCC=-0.000165, ΔSCI=-0.001994, ΔSPS=-0.000453` だった。

全85のうち pilot12 を分離して再集計した（case-level の記述的分割であり、
pilot と nonpilot の family は共有されるため family-held-out 評価ではない）。
nonpilot73 の baseline500 対差は次の通り。

| 条件 | ΔMCC | ΔSCI | ΔSPS | 緩い0.005基準 |
|---|---:|---:|---:|---|
| LPC projected 200 | −0.002504 | −0.005777 | −0.000806 | SCI が不通過 |
| LPV projected 300 | −0.001460 | −0.001661 | −0.000660 | 3指標とも範囲内 |

これは case-level の記述的分割であり、新しい bootstrap や推論を実施したもの
ではない。

baseline500 対の品質 worst individual は、LPC projected 200 で MCC
`−0.049483`（RF00017-2）、SCI `−0.234212`（RF00017-4）、SPS `−0.036216`
（RF00031-1）、LPV projected 300 で MCC `−0.044740`（RF00059-0）、SCI
`−0.255100`（RF00017-3）、SPS `−0.016584`（RF00436-1）だった。平均だけで
等価と判断できない理由である。

## 理論上の位置づけ

`column-bound` は UB/stopping を変えるだけで、stop までの merge trajectory は
同じであるため、5倍の品質改善を期待する対象ではない。`tight-alignment` は
certified oracle と追加の feasible repair を変える歴史的比較対象だが、
O(M log L) のため、O(L+M) の strict-linear 提案 default として扱わない。
BeamAlign と LinearNuss は変更していない。

なお、これは production の global convergence や lower cap の証明ではない。
表の `gap`（`dd_final_gap_last`）は merge ごとの平均ではなく、最後の
progressive merge における certified gap `bestUB−LB` である。
`dd_summary.violated` は最後の内部 iterate の dual disagreement であり、最終
output の構造が invalid だという意味ではない。DD が選択する best-feasible merge
solution と、その後に別途 decode される final consensus structure を区別し、
内部 violation、certified violation、最終構造の feasibility を別々に読む。

bulk の `wall`/`DDwall` は run.py の jobs=16、metrics-jsonl 有効な throughput
観測であり、固定CPUでの純粋な attribution ではない。metrics overhead を外した
比較は、同じケースを interleave した serial metrics-off arm によって別途記録する。

## 成果物・検証コマンド

主スクリプトは [benchmark_dd_budget.py](../benchmarks/benchmark_dd_budget.py)、
テストは [test_benchmark_dd_budget.py](../benchmarks/test_benchmark_dd_budget.py)。
生成 config と run manifest は次の root に保存した。

`benchmarks/results-murlet-dataset1-dd-budget-arm64-v272/`

主要な監査ファイルは `plan.json`、`dd_budget_report.json`、
`dd_budget_per_condition.csv`、`dd_budget_pairs.csv`、各 scope の `manifest.json`
である。future-generated config からは旧 `comparison`/`luna_relaxed` metadata を
除去し、既に実行済みの pilot config は上書きしていない。

```text
python3 -m unittest discover -s benchmarks -p 'test_benchmark_dd_budget.py' -v
python3 benchmarks/benchmark_dd_budget.py summarize \
  --results-root benchmarks/results-murlet-dataset1-dd-budget-arm64-v272
cmake --build build-arm64 -j12
ctest --test-dir build-arm64 --output-on-failure
sha256sum build-arm64/src/dafs
find benchmarks/results-murlet-dataset1-dd-budget-arm64-v272 -name manifest.json -print | sort
sha256sum benchmarks/configs/dd-budget-*.json \
  benchmarks/configs/luna-relaxed.json benchmarks/data/murlet/manifest.json
python3 benchmarks/benchmark_dd_budget.py serial --metrics-off \
  --results-root benchmarks/results-murlet-dataset1-dd-budget-arm64-v272/serial-cpu19-metrics-off \
  --cpu 19 --repetitions 5
```

build/ctest は現行 dirty source tree で実行し、10/10 test pass だった。production
source、default、BeamAlign/LinearNuss をこの評価のためには変更していない。

## metrics-off 固定CPU timing

品質選択後に性能だけを見てケースを選ばないよう、入力 alignment length で固定した
short/long の `RF00169-2` と `RF00017-1` を、pilot 後・timing 前に選んだ。これは
厳密な pilot 前登録でも性能ベースの選択でもない。CPU19、OMP=1、5 repetitions、
arm 順を interleave し、`--metrics-jsonl` を付けない別 root
`serial-cpu19-metrics-off/` で80 runを実行した。全80 prediction SHA256 が bulk の
対応 arm（rep-01）と一致し、status failure は0件だった。
先行して完了していた `serial-cpu19/` は `metrics-on` の補足 throughput 観測として
別ラベルで保存し、metrics-off の attribution と混ぜていない。

下表の speedup は reference wall / candidate wall。各セルは5回の per-case median、
括弧内は2ケース×5回の ratio の mean / median である。秒の aggregate を1ケースの
時間とは解釈しない。metrics-off なので DDwall は測定せず、品質表の DDwall は
jobs=16・metrics-on bulk の throughput 観測としてのみ使う。

| 条件・candidate | reference | RF00017-1 | RF00169-2 | aggregate speedup (mean / median) |
|---|---|---:|---:|---:|
| LPC×LinearCONTRAlign projected 200 | baseline 500 | 1.606 / 1.607 CPU | 1.788 / 1.798 CPU | wall 1.698 / 1.684; CPU 1.702 / 1.705 |
| LPC×LinearCONTRAlign projected 200 | baseline 200 | 1.001 / 1.001 CPU | 0.988 / 0.988 CPU | wall 0.997 / 1.000; CPU 0.996 / 0.993 |
| LPC×LinearCONTRAlign projected 200 | baseline 100 | 0.794 / 0.793 CPU | 0.726 / 0.726 CPU | wall 0.758 / 0.761; CPU 0.758 / 0.759 |
| LPV×LinearProbCons projected 300 | baseline 500 | 1.383 / 1.384 CPU | 1.447 / 1.450 CPU | wall 1.414 / 1.410; CPU 1.415 / 1.410 |
| LPV×LinearProbCons projected 300 | baseline 300 | 1.012 / 1.012 CPU | 0.986 / 0.993 CPU | wall 1.001 / 0.999; CPU 1.002 / 1.000 |
| LPV×LinearProbCons projected 300 | baseline 100 | 0.613 / 0.614 CPU | 0.514 / 0.514 CPU | wall 0.564 / 0.567; CPU 0.564 / 0.567 |

したがってこの固定した2ケースでは projected arm は baseline500 より速いが、
same-budget baseline とほぼ同速であり、baseline100 より速いわけではない。この
timing 結果を品質結果と合わせても、maxiter=100≈500 や普遍的な5x改善は主張しない。
