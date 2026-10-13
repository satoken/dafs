# Murlet dataset1 条件スイープ結果

Murlet `dataset1` 全85ケースについて、前回の条件探索範囲をarm64版DAFSで再評価した。8構成（VP/CONTRAfoldのAlifoldあり・なしを含む）、320条件、合計27,200 runである。各条件は各ケース1回、平均値は85ケースのmacro averageである。SCIはViennaRNA 2.7.2で計算した。

条件の順位は、平均MCCを第一キー、平均SCIを第二キー、平均SPSを第三キーとする辞書式順序で決めた。MCCに許容差は設けていない。全27,200 runで `status=success`、`metrics_valid=True`、`quality_complete=True` だった。

## 組み合わせ別の最良条件

| 組み合わせ | 選択条件 | phase | DD RIBOSUM | 最終RIBOSUM | SPS | alignment PPV | 構造感度 | 構造PPV | MCC | CBP F1 | SCI |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ViennaRNA × ProbCons | `vp_alifold_w2` | parameter | 0.02 | 0.0 | 0.7657 | 0.7714 | 0.6559 | 0.7320 | **0.6812** | 0.5790 | 0.6622 |
| CONTRAfold × CONTRAlign | `cc_alifold_fth03` | parameter | 0.1 | 0.0 | 0.7732 | 0.7810 | 0.6714 | 0.7537 | **0.6999** | 0.5972 | 0.6805 |
| lpc × LinearAlign | `cl_linear_ath02` | parameter | 0.5 | 0.01 | 0.7572 | 0.7752 | 0.6711 | 0.7156 | **0.6851** | 0.5749 | 0.6638 |
| lpv × LinearProbCons | `lpv_probcons_fpct05` | parameter | 0.3 | 0.0 | 0.7564 | 0.7667 | 0.6605 | 0.6858 | **0.6635** | 0.5543 | 0.6484 |
| lpc × LinearCONTRAlign | `lpc_contralign_ath02` | parameter | 0.3 | 0.005 | 0.7530 | 0.7658 | 0.6587 | 0.7250 | **0.6833** | 0.5754 | 0.6513 |
| lpv × LinearAlign | `lpv_linear_fpct05` | parameter | 0.3 | 0.0 | 0.7751 | 0.7855 | 0.6690 | 0.6905 | **0.6702** | 0.5710 | 0.6630 |

## phase別のMCC最良条件

| 組み合わせ | base | parameter | DD RIBOSUM | 最終RIBOSUM |
|---|---|---|---|---|
| ViennaRNA × ProbCons | `vp_alifold_base` (MCC 0.6746, SCI 0.6819, SPS 0.7559) | `vp_alifold_w2` (MCC 0.6812, SCI 0.6622, SPS 0.7657) | `vp_alifold_dd_0p075` (MCC 0.6771, SCI 0.6810, SPS 0.7542) | `vp_alifold_final_0p01` (MCC 0.6756, SCI 0.6819, SPS 0.7559) |
| CONTRAfold × CONTRAlign | `cc_alifold_base` (MCC 0.6921, SCI 0.7036, SPS 0.7688) | `cc_alifold_fth03` (MCC 0.6999, SCI 0.6805, SPS 0.7732) | `cc_alifold_dd_0` (MCC 0.6936, SCI 0.7097, SPS 0.7741) | `cc_alifold_final_0p005` (MCC 0.6924, SCI 0.7036, SPS 0.7688) |
| lpc × LinearAlign | `cl_linear_base` (MCC 0.6803, SCI 0.6576, SPS 0.7564) | `cl_linear_ath02` (MCC 0.6851, SCI 0.6638, SPS 0.7572) | `cl_linear_dd_0p3` (MCC 0.6843, SCI 0.6694, SPS 0.7746) | `cl_linear_final_0p01` (MCC 0.6803, SCI 0.6576, SPS 0.7564) |
| lpv × LinearProbCons | `lpv_probcons_base` (MCC 0.6613, SCI 0.6686, SPS 0.7537) | `lpv_probcons_fpct05` (MCC 0.6635, SCI 0.6484, SPS 0.7564) | `lpv_probcons_dd_0p1` (MCC 0.6626, SCI 0.6641, SPS 0.7652) | `lpv_probcons_final_0p005` (MCC 0.6573, SCI 0.6686, SPS 0.7537) |
| lpc × LinearCONTRAlign | `lpc_contralign_base` (MCC 0.6771, SCI 0.6436, SPS 0.7542) | `lpc_contralign_ath02` (MCC 0.6833, SCI 0.6513, SPS 0.7530) | `lpc_contralign_dd_0p5` (MCC 0.6822, SCI 0.6367, SPS 0.7457) | `lpc_contralign_final_0p01` (MCC 0.6790, SCI 0.6436, SPS 0.7542) |
| lpv × LinearAlign | `lpv_linear_base` (MCC 0.6629, SCI 0.6757, SPS 0.7654) | `lpv_linear_fpct05` (MCC 0.6702, SCI 0.6630, SPS 0.7751) | `lpv_linear_dd_0p15` (MCC 0.6662, SCI 0.6779, SPS 0.7808) | `lpv_linear_final_0p005` (MCC 0.6585, SCI 0.6757, SPS 0.7654) |

## 選択条件の全引数

### ViennaRNA × ProbCons — `vp_alifold_w2`

```text
-a ProbCons -s Vienna --alifold --dynamic-cbp --ribosum-weight=0.02 --final-ribosum-weight=0 --fold-th=0.3 --max-iter=80 --weight=2
```

### CONTRAfold × CONTRAlign — `cc_alifold_fth03`

```text
-a CONTRAlign -s CONTRAfold --alifold --dynamic-cbp --ribosum-weight=0.1 --final-ribosum-weight=0 --weight=2 --max-iter=80 --fold-th=0.3
```

### lpc × LinearAlign — `cl_linear_ath02`

```text
-a LinearAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.5 --final-ribosum-weight=0.01 --max-iter=80 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100 --align-th=0.02
```

### lpv × LinearProbCons — `lpv_probcons_fpct05`

```text
-a LinearAlign-ProbConsRNA -s lpv --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --weight=4 --max-iter=80 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100 --fold-pct=0.5
```

### lpc × LinearCONTRAlign — `lpc_contralign_ath02`

```text
-a LinearAlign-CONTRAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0.005 --align-pct=0.5 --max-iter=80 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100 --align-th=0.02
```

### lpv × LinearAlign — `lpv_linear_fpct05`

```text
-a LinearAlign -s lpv --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --align-pct=0.5 --max-iter=80 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100 --fold-pct=0.5
```

## LinearTurboFold baseline

LinearTurboFoldは条件スイープ対象ではなく、前回と同じ固定条件（85ケース、hmm/CKY beam 100、3反復、threshold 0.3）をbaselineとして併記する。

| 方法 | SPS | alignment PPV | 構造感度 | 構造PPV | MCC | CBP F1 | SCI |
|---|---:|---:|---:|---:|---:|---:|---:|
| LinearTurboFold | 0.8130 | 0.8075 | 0.6757 | 0.5771 | 0.6134 | 0.5381 | 0.5301 |

## 実験ファイル

- [スイープ設定](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/configs/murlet-dataset1-condition-sweep-arm64-v272.json)
- [集計CSV](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/results-murlet-dataset1-condition-sweep-arm64-v272/summary.csv)
- [LinearTurboFold結果](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/results-linearturbofold-murlet-dataset1-arm64-v272/summary.csv)
