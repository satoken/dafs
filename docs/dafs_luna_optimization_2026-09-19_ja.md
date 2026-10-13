# DAFS Luna最適化ベンチマーク・比較レポート（2026-09-19）

> decoder-only と DD/tight-alignment の full 680 比較は完了。tight-alignment は opt-in の実験結果として記録する。

## 目的と実験範囲

- 入力設定は [`murlet-dataset1-buffer-reuse-arm64-v272.json`](../benchmarks/configs/murlet-dataset1-buffer-reuse-arm64-v272.json) から生成する。
- decoder-only（既定）は、新しいDDフラグを付けない。4条件 × `max-iter=100,1000` × Murlet dataset1 85例 = 680 run。
- baseline binary は `/tmp/dafs-before-luna-20260919`、optimized binary は `build-arm64/src/dafs`。
- decoder 比較に使用した binary SHA-256 は、baseline `7c6cfe1d01fe5a37404dfb7f4193fc018c056369d6170f6289133dfa0e5151c3`、frozen optimized `/tmp/dafs-luna-decoder-20260919` `49f9abb54daef26701ef0b396192b113c47ee851ddf6c86d23f5d687eee43432`。
- 既定の比較は、直前の buffer-reuse 結果を historical baseline として、optimized の fresh result tree と同一の condition/dataset/rep を対応付ける。
- tight-alignment は `--dd-tight-alignment` で明示的に有効化する。同じ 680 key を比較するが、prediction・quality・DD path の変更は許容し、欠落 run と実行/metrics error は許容しない。
- full run は `--jobs 16`。表の時間は DAFS process の計測値で、外部 SCI 評価時間を含まない。

## 再現コマンド

生成は設定ファイルだけを書き、build と benchmark execution は行わない。

```sh
python3 benchmarks/benchmark_luna_optimization.py generate \
  --config-dir benchmarks/configs \
  --results-root benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272 \
  --baseline-results /tmp/dafs-luna-decoder-baseline \
  --optimized-results benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272 \
  --baseline-binary /tmp/dafs-before-luna-20260919 \
  --optimized-binary build-arm64/src/dafs

python3 benchmarks/run.py --config benchmarks/configs/luna-decoder-optimized.json \
  --jobs 16 --binary /tmp/dafs-luna-decoder-20260919

python3 benchmarks/benchmark_luna_optimization.py compare \
  --config benchmarks/configs/luna-decoder-optimized.json \
  --baseline-results benchmarks/results-murlet-dataset1-buffer-reuse-arm64-v272 \
  --optimized-results benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272 \
  --baseline-kind historical \
  --compare-traces \
  --report-dir benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272/comparison
```

Tight-alignment branch:

```sh
python3 benchmarks/benchmark_luna_optimization.py generate \
  --dd-tight-alignment \
  --config-dir benchmarks/configs \
  --results-root benchmarks/results-murlet-dataset1-luna-tight-arm64-v272 \
  --baseline-results /tmp/dafs-luna-tight-baseline \
  --optimized-results benchmarks/results-murlet-dataset1-luna-tight-arm64-v272

python3 benchmarks/run.py --config benchmarks/configs/luna-tight-optimized.json \
  --jobs 16 --binary /tmp/dafs-luna-decoder-20260919
python3 benchmarks/benchmark_luna_optimization.py compare \
  --config benchmarks/configs/luna-tight-optimized.json \
  --baseline-results benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272 \
  --optimized-results benchmarks/results-murlet-dataset1-luna-tight-arm64-v272 \
  --baseline-kind fresh \
  --report-dir benchmarks/results-murlet-dataset1-luna-tight-arm64-v272/comparison
```

Fresh timing pilot（例）では `--dataset-suffix=-0 --limit 3` を両 binary に同じく適用し、main が matched case を baseline/new 交互に実行する。driver は key matching と集計を行うが、熱・負荷・scheduler contention の制御を仮定しない。

## 出力と判定

各 compare report は次を出力する。

- `luna_optimization_per_condition.csv`: condition ごとの paired count、prediction byte equality、quality dict equality、run-summary non-time equality、MCC/SCI/SPS、DD iterations、stop reason、violations、final/certified gap、DD/wall/CPU/RSS の平均。
- `luna_optimization_report.json`: per-run の対応付け、欠落/余剰 run、errors、equality counts、per-condition 集計、provenance、制限事項。

速度比は `baseline / optimized`、時間削減率は `100 * (baseline - optimized) / baseline` とする。RSS は `100 * (optimized - baseline) / baseline` の変化率であり、速度比とは解釈しない。

`dd_seconds` は DAFS の `steady_clock` による DD wall time であり CPU time ではない。したがって、DD time を CPU-only の改善、または contention 非依存の改善とは主張しない。historical baseline と fresh optimized の比較は paired key であっても同一時刻の制御実験ではない。

decoder-only では prediction byte、quality dict、run-summary の time fields を除く内容の不一致（`--compare-traces` 指定時は正規化後の全event不一致）で終了コード non-zero とする。tight-alignment ではそれらの出力差を報告するが、欠落 run・不成功 result・metrics parse/validation error は non-zero とする。

## 理論上の位置付け

参照資料は Stanford EE364b の [Subgradient Methods notes](https://stanford.edu/class/ee364b/lectures/subgrad_method_notes.pdf)（§4.2 Polyak step size choice with estimated (f^\star)、§6.2 projected subgradient for the dual problem）。今回の解釈は次の通り。

- exact subgradient と、beam による近似 decoder が返す方向は同一視しない。
- feasible primal の lower bound は、dual optimum と等しいとは限らない。したがって feasible repair や certified track の値は、観測された bound として報告する。
- constant Polyak step を採用したことだけから、一般の収束保証を主張しない。実験では stop reason、violations、final gap、certified gap を併記する。

DD tight-alignment の Fenwick 実装は (O(m\log L)) の明示的な非 strict-linear 実装であり、既定の decoder-only path には影響させない。既存の exact pricing と tolerance は維持する、という前提で比較する。

## 実装概要と保証範囲

- BeamAlign は sorted compact reusable layers を使い、`nth_element` 後に original key order と tie order を復元する。既定の arithmetic と tie semantics は変更しない。
- LinearNussinov は generation stamps、suffix query sweep、root cache を使って再利用する。こちらも既定の arithmetic と tie semantics は変更しない。
- 「length-linear」は beam=100 固定の今回の実験範囲に限る。vector insertion は layer ごとの worst case が `O(b²)` であり、任意の beam に対する普遍的な速度保証ではない。

実装検証は ctest `9/9` pass、standalone Nussinov の ASAN+UBSAN pass（`ASAN_OPTIONS=detect_leaks=0`）である。環境では LSan が ptrace により利用できないため、leak sanitizer の結果は主張しない。

## 実装局所ベンチマーク（提供値）

BeamAlign の reused map → compact vector の ARM64/O3 3-run median 測定値。速度比はここでも old/new とする。

| shape | beam | calls | old ms | new ms | speedup |
|---|---:|---:|---:|---:|---:|
| 60×60 | 100 | 100 | 22.62 | 10.20 | 2.22× |
| 100×100 | 100 | 30 | 19.96 | 9.95 | 2.01× |
| 100×100 | 1000 | 10 | 6.78 | 3.60 | 1.88× |
| 400×400 | 100 | 5 | 26.37 | 13.81 | 1.91× |
| 400×400 | 1000 | 5 | 58.21 | 41.49 | 1.40× |
| 250×250 | 0 | 10 | 44.53 | 28.03 | 1.59× |
| 400×400 | 0 | 5 | 59.64 | 41.61 | 1.43× |

この表の範囲は 1.40–2.22×。全体 benchmark の結論とは分けて扱う。

## Controlled paired timing anchor

`RF00169-2 / cl_linear_ath02_iter1000` を、full decoder run 完了後に idle を挟んで交互実行した metrics-off A/B 3 pairs。binary は上記の baseline/frozen optimized SHA-256 に固定した。

- before wall: `7.157363, 7.140274, 7.116635` s、CPU: `7.146083, 7.137689, 7.113949` s
- after wall: `4.642156, 4.608320, 4.578695` s、CPU: `4.638606, 4.600916, 4.577327` s
- median wall: `7.140274→4.608320` s（約 `1.55×`）、median CPU: `7.137689→4.600916` s（約 `1.55×`）
- timing は DAFS process のみで、外部 SCI 評価を含まない。prediction identity は full comparison で確認する。

## Full 680 結果

### decoder-only

Full 680 は完走した。`expected=680, paired=680, errors=0, missing=0, unexpected=0`。prediction byte、quality dict、run-summary の非時間部分、`--compare-traces` の正規化後 event がすべて `680/680` 一致した。

速度・RSS の平均は 85 dataset cases。MCC/SCI/SPS は baseline と optimized が同値なので optimized 値を記載する。`wall`、`DD wall`、`CPU` は old→new、RSS は optimized の変化率。

| condition | equal (pred/quality/summary/trace) | MCC | SCI | SPS | wall s | wall × | DD wall s | DD × | CPU s | RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `cl_linear_ath02_iter100` | 85/85/85/85 | .6852 | .6792 | .7553 | 4.4362→2.9699 | 1.494 | 2.6007→1.2808 | 2.031 | 4.4199→2.9542 | −2.33% |
| `cl_linear_ath02_iter1000` | 85/85/85/85 | .6913 | .7445 | .7559 | 22.7761→12.5348 | 1.817 | 20.8625→10.7234 | 1.946 | 22.7155→12.4722 | −2.14% |
| `lpc_contralign_ath02_iter100` | 85/85/85/85 | .6857 | .6623 | .7536 | 4.3576→2.9311 | 1.487 | 2.4523→1.2223 | 2.006 | 4.3424→2.9128 | −1.81% |
| `lpc_contralign_ath02_iter1000` | 85/85/85/85 | .6933 | .7184 | .7529 | 21.4292→11.6051 | 1.847 | 19.4637→9.7788 | 1.990 | 21.3702→11.5446 | −1.92% |
| `lpv_linear_fpct05_iter100` | 85/85/85/85 | .6694 | .6809 | .7728 | 4.6162→3.3025 | 1.398 | 2.6808→1.4820 | 1.809 | 4.6007→3.2875 | −2.62% |
| `lpv_linear_fpct05_iter1000` | 85/85/85/85 | .6747 | .7642 | .7686 | 24.2474→14.2131 | 1.706 | 22.3071→12.4067 | 1.798 | 24.1907→14.1494 | −2.46% |
| `lpv_probcons_fpct05_iter100` | 85/85/85/85 | .6635 | .6609 | .7573 | 4.4011→3.1306 | 1.406 | 2.5107→1.3766 | 1.824 | 4.3848→3.1147 | −2.52% |
| `lpv_probcons_fpct05_iter1000` | 85/85/85/85 | .6703 | .7258 | .7566 | 23.1158→13.1225 | 1.762 | 21.1982→11.2958 | 1.877 | 23.0535→13.0618 | −2.58% |

`DD ×` は baseline DD wall / optimized DD wall、すなわち速度比。`dd_seconds` は `steady_clock` wall time で CPU time ではない。

DD diagnostics are identical baseline→optimized for every condition; the following are the common per-case means/counts from the report.

| condition | DD iterations | stop reasons (`n`) | violations/case | certified violations/case | final gap mean | certified gap mean |
|---|---:|---|---:|---:|---:|---:|
| `cl_linear_ath02_iter100` | 893.6 | beam_stationary 16; max_iterations 749 | 280.8 | 410.6 | 13.9187 | 13.9187 |
| `cl_linear_ath02_iter1000` | 8264.9 | beam_stationary 71; certified_gap 19; max_iterations 675 | 107.0 | 275.5 | 10.1917 | 10.1917 |
| `lpc_contralign_ath02_iter100` | 894.8 | beam_stationary 17; max_iterations 748 | 225.4 | 280.3 | 8.6857 | 8.6857 |
| `lpc_contralign_ath02_iter1000` | 8061.3 | beam_stationary 90; certified_gap 25; max_iterations 650 | 86.8 | 174.3 | 6.3919 | 6.3919 |
| `lpv_linear_fpct05_iter100` | 899.3 | beam_stationary 3; max_iterations 762 | 390.5 | 456.6 | 11.0096 | 11.0096 |
| `lpv_linear_fpct05_iter1000` | 8402.6 | beam_stationary 73; certified_gap 3; max_iterations 689 | 118.3 | 267.7 | 7.7483 | 7.7483 |
| `lpv_probcons_fpct05_iter100` | 898.9 | beam_stationary 9; max_iterations 756 | 318.6 | 354.3 | 9.0329 | 9.0329 |
| `lpv_probcons_fpct05_iter1000` | 8252.0 | beam_stationary 92; certified_gap 4; max_iterations 669 | 95.6 | 220.1 | 6.9075 | 6.9075 |

Exact source artifacts: [`luna_optimization_report.json`](../benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272/comparison/luna_optimization_report.json) and [`luna_optimization_per_condition.csv`](../benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272/comparison/luna_optimization_per_condition.csv).

### tight-alignment

Full DD 680 は全 run が成功・valid・quality complete で完了した。`expected=680, paired=680, errors=0, missing=0, unexpected=0`。`--compare-traces` は tight branch では指定していない（`compare_traces=false`）。tight は出力変更を許容する opt-in 比較なので、prediction byte は `304/680` が一致、`376/680` が変更、quality dict は `306/680` が一致、summary の非時間部分は `3/680` が一致した。したがって、これは decoder-only の回帰判定ではなく、DD path の出力変更を含む実験結果である。

以下は 85 dataset cases の平均。tight 比較の `old→new` は、同じ frozen optimized decoder binary における DD-off → 同じ binary の DD-on（tight-alignment）であり、prior buffer-reuse/historical baseline ではない。`DD iter` は 1 case の 9 merge 合計、`wall ×` は old/new、`wall Δ` は `100*(new-old)/old` である。

| condition | equal (pred/quality/summary) | MCC | SCI | SPS | DD iter | DD wall s | wall s | wall × | wall Δ | CPU s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `cl_linear_ath02_iter100` | 24/24/1 | .6852→.6839 | .6792→.6878 | .7553→.7560 | 893.6→893.1 | 1.281→1.627 | 2.970→3.381 | .878 | +13.8488% | 2.954→3.364 |
| `cl_linear_ath02_iter1000` | 57/57/0 | .6913→.6907 | .7445→.7441 | .7559→.7560 | 8264.9→7716.5 | 10.723→12.681 | 12.535→14.513 | .864 | +15.7782% | 12.472→14.438 |
| `lpc_contralign_ath02_iter100` | 28/29/1 | .6857→.6895 | .6623→.6689 | .7536→.7557 | 894.8→895.0 | 1.222→1.510 | 2.931→3.327 | .881 | +13.5225% | 2.913→3.310 |
| `lpc_contralign_ath02_iter1000` | 55/55/0 | .6933→.6931 | .7184→.7208 | .7529→.7525 | 8061.3→7485.7 | 9.779→11.733 | 11.605→13.574 | .855 | +16.9683% | 11.545→13.503 |
| `lpv_linear_fpct05_iter100` | 25/25/1 | .6694→.6676 | .6809→.6920 | .7728→.7725 | 899.3→899.0 | 1.482→1.817 | 3.302→3.698 | .893 | +11.9803% | 3.288→3.679 |
| `lpv_linear_fpct05_iter1000` | 47/48/0 | .6747→.6734 | .7642→.7679 | .7686→.7687 | 8402.6→7966.3 | 12.407→15.225 | 14.213→17.133 | .830 | +20.5461% | 14.149→17.066 |
| `lpv_probcons_fpct05_iter100` | 21/21/0 | .6635→.6642 | .6609→.6736 | .7573→.7568 | 898.9→898.4 | 1.377→1.666 | 3.131→3.448 | .908 | +10.1428% | 3.115→3.429 |
| `lpv_probcons_fpct05_iter1000` | 47/47/0 | .6703→.6709 | .7258→.7302 | .7566→.7562 | 8252.0→7635.6 | 11.296→12.877 | 13.122→14.657 | .895 | +11.6953% | 13.062→14.587 |

DD stop reasons are counted across 765 merges per condition (`85×9`). `viol/case` and `cert-viol/case` are the corresponding report means; gaps are the report mean final gap.

| condition | stop reasons old→new (`n`) | violations/case old→new | certified violations/case old→new | final gap old→new | certified gap old→new |
|---|---|---:|---:|---:|---:|
| `cl_linear_ath02_iter100` | beam_stationary 16; max_iterations 749 → beam_stationary 18; certified_gap 3; max_iterations 744 | 280.8→278.2 | 410.6→323.9 | 13.9187→4.7368 | 13.9187→4.7368 |
| `cl_linear_ath02_iter1000` | beam_stationary 71; certified_gap 19; max_iterations 675 → beam_stationary 86; certified_gap 122; max_iterations 557 | 107.0→107.2 | 275.5→147.8 | 10.1917→.4074 | 10.1917→.4074 |
| `lpc_contralign_ath02_iter100` | beam_stationary 17; max_iterations 748 → beam_stationary 13; certified_gap 4; max_iterations 748 | 225.4→233.9 | 280.3→263.0 | 8.6857→3.1781 | 8.6857→3.1781 |
| `lpc_contralign_ath02_iter1000` | beam_stationary 90; certified_gap 25; max_iterations 650 → beam_stationary 100; certified_gap 138; max_iterations 527 | 86.8→87.7 | 174.3→91.8 | 6.3919→.3000 | 6.3919→.3000 |
| `lpv_linear_fpct05_iter100` | beam_stationary 3; max_iterations 762 → beam_stationary 5; max_iterations 760 | 390.5→380.0 | 456.6→451.0 | 11.0096→4.6319 | 11.0096→4.6319 |
| `lpv_linear_fpct05_iter1000` | beam_stationary 73; certified_gap 3; max_iterations 689 → beam_stationary 98; certified_gap 72; max_iterations 595 | 118.3→115.6 | 267.7→141.8 | 7.7483→.4074 | 7.7483→.4074 |
| `lpv_probcons_fpct05_iter100` | beam_stationary 9; max_iterations 756 → beam_stationary 8; max_iterations 757 | 318.6→315.4 | 354.3→336.9 | 9.0329→3.2442 | 9.0329→3.2442 |
| `lpv_probcons_fpct05_iter1000` | beam_stationary 92; certified_gap 4; max_iterations 669 → beam_stationary 128; certified_gap 99; max_iterations 538 | 95.6→96.3 | 220.1→122.1 | 6.9075→.2980 | 6.9075→.2980 |

max-iter 1000 では DD iteration 数が cl `6.6357%`、lpc-contra `7.1410%`、lpv-linear `5.1921%`、lpv-prob `7.4699%` 減少した。一方、full wall はそれぞれ `+15.7782%`, `+16.9683%`, `+20.5461%`, `+11.6953%` で、全体 runtime は回帰した。max-iter 100 は iteration 差がほぼ無く、lpc-contra のみ僅かな増加（894.8→895.0）だった。`certified_gap` stop は max-iter 1000 の4条件で `19→122`, `25→138`, `3→72`, `4→99`（各 765 merge）となり、gap 縮小の機構を示すが、時間回帰を相殺しない。

max-iter 1000 の `certified_alignment_repair_seconds`（9 merge 合計）の per-case 平均は、cl `1.1864 s`（DD の `9.36%`）、lpc-contra `1.0187 s`（`8.68%`）、lpv-linear `1.3127 s`（`8.62%`）、lpv-prob `1.2534 s`（`9.73%`）。repair は overhead の一部だが、全 overhead を repair だけには帰属しない（Fenwick oracle 2 系統と support-union も追加されている）。将来の候補として、既存の validUB/exact pricing を維持したまま feasible repair の next candidate を安全に amortize/cache できる余地はあるが、本報告では実装したとは主張しない。

追加の serial metrics-off A/B（`RF00169-2 / cl_linear_ath02_iter1000`、DD-off/on を交互に3 pairs）は、off wall `4.641752 s` → on `5.537914 s`（約 `+19.31%`）、off CPU `4.630351 s` → on `5.536750 s`（約 `+19.58%`）の median だった。各群内は安定していたが、`dd_seconds` は `steady_clock` wall time なので CPU-only、または contention 非依存の改善とは解釈しない。出力変更は opt-in 比較として許容した。

結論として、tight-alignment は gap と certified stop を改善する optional proof-preserving experiment path として保持するが、total runtime が全条件で増えたため、既定 ON は推奨しない。

Exact source artifacts: [`luna_optimization_report.json`](../benchmarks/results-murlet-dataset1-luna-tight-arm64-v272/comparison/luna_optimization_report.json) and [`luna_optimization_per_condition.csv`](../benchmarks/results-murlet-dataset1-luna-tight-arm64-v272/comparison/luna_optimization_per_condition.csv).

## 再現 artifact

- Generated configs: `benchmarks/configs/luna-decoder-*.json`, `benchmarks/configs/luna-tight-*.json`
- Decoder CSV/JSON: `benchmarks/results-murlet-dataset1-luna-decoder-arm64-v272/comparison/`
- Tight CSV/JSON: `benchmarks/results-murlet-dataset1-luna-tight-arm64-v272/comparison/`
