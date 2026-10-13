# Murlet dataset1 全85ケース・ベンチマーク

各組み合わせについて、Murlet dataset1 の13ケースでMCCが最良だった条件を固定し、dataset1全85ケースを1回ずつ実行した。DAFSは6条件、LinearTurboFoldは1条件で、合計595 runを評価した。平均値はケースごとの単純平均（macro average）である。SCIはarm64版ViennaRNA 2.7.2で計算した。

## 精度

| 組み合わせ | run | SPS | alignment PPV | 構造 sensitivity | 構造 PPV | MCC | CBP F1 | SCI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CONTRAfold × CONTRAlign | 85 | 0.7649 | 0.7755 | 0.6754 | 0.7060 | **0.6828** | 0.5764 | 0.6867 |
| lpc × LinearAlign | 85 | 0.7564 | 0.7748 | 0.6683 | 0.7083 | 0.6803 | 0.5699 | 0.6576 |
| lpc × LinearCONTRAlign | 85 | 0.7542 | 0.7648 | 0.6530 | 0.7191 | 0.6771 | 0.5723 | 0.6436 |
| ViennaRNA × ProbCons | 85 | 0.7559 | 0.7690 | 0.6569 | 0.7166 | 0.6746 | 0.5643 | 0.6819 |
| lpv × LinearAlign | 85 | 0.7654 | 0.7805 | 0.6806 | 0.6645 | 0.6629 | 0.5582 | 0.6757 |
| lpv × LinearProbCons | 85 | 0.7537 | 0.7673 | 0.6763 | 0.6660 | 0.6613 | 0.5511 | 0.6686 |
| LinearTurboFold | 85 | 0.8130 | 0.8075 | 0.6757 | 0.5771 | 0.6134 | 0.5381 | 0.5301 |

全595 runで、DAFSは `status=success`、`metrics_valid=True`、`quality_complete=True` だった。LinearTurboFoldも85/85件が成功し、全ケースで参照との精度を計算できた。

## 固定した条件

DAFSの共通条件は `--max-iter=80`、`OMP_NUM_THREADS=1`。線形条件では `--align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100` とした。全引数は [DAFS設定](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/configs/murlet-dataset1-best-arm64-v272.json) に保存している。

| 条件ID | 主な追加条件 |
|---|---|
| `vp_alifold_mcc` | ProbCons × Vienna、RNAalifold有効、DD RIBOSUM 0.02、`fold-th=0.3` |
| `cc_noalifold_mcc` | CONTRAlign × CONTRAfold、DD RIBOSUM 0.1、`fold-pct=0` |
| `cl_linear_mcc` | LinearAlign × lpc、DD RIBOSUM 0.5、最終RIBOSUM 0.01、`align-th=0.005` |
| `lpv_probcons_base` | LinearAlign-ProbConsRNA × lpv、DD RIBOSUM 0.3、`weight=4` |
| `lpc_contralign_alignpct05` | LinearAlign-CONTRAlign × lpc、DD RIBOSUM 0.3、最終RIBOSUM 0.005、`align-pct=0.5` |
| `lpv_linear_mcc_apct05` | LinearAlign × lpv、DD RIBOSUM 0.3、`align-pct=0.5` |

LinearTurboFoldは `hmm_beam=100`、`cky_beam=100`、3反復、`min_helix_length=3`、`pk_iterations=1`、`threshold=0.3` とした。設定は [LinearTurboFold設定](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/configs/linearturbofold-murlet-dataset1-arm64-v272.json) に保存している。

## 結果ファイル

- [DAFS summary.csv](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/results-murlet-dataset1-best-arm64-v272/summary.csv)
- [LinearTurboFold summary.csv](/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/results-linearturbofold-murlet-dataset1-arm64-v272/summary.csv)

13ケースだけで選んだ条件をdataset1全体へ適用すると、MCCはnonlinearのCONTRAfold × CONTRAlignが最良だった。LinearTurboFoldはSPSとalignment PPVは最良だが、構造PPV、MCC、SCIは低く、アラインメント精度と共通構造精度のトレードオフが全85ケースでも確認された。
