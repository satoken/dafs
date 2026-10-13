# Murlet dataset1 LinearAlifold 条件スイープ

85ケース、480条件、40800 run。
全条件で `--max-iter=100`、linear系の5種類のbeamを100に固定した。
順位は85ケースのmacro平均MCC、SCI、SPSの辞書式順序で決めた。
全runで `status=success`、`metrics_valid=True`、`quality_complete=True` を確認した。
全runの開始時刻と実行ファイルのパスをmanifestと照合し、バイナリと設定のSHA256も確認した。
評価バイナリ SHA256: `b5e3542853b6e368a88d483787365745d1805adedf36ebcb2a74e8ddbbd6fee7`。
設定 SHA256: `a4243213dff3595a81c2f01166c8512cd0eb1430f3c891f39f8ad74d34fb2efa`。実行開始: `2026-09-19T09:12:24.209757+00:00`。
同じケースを条件選択と結果報告に使う探索結果であり、独立したテスト集合の精度ではない。

## 構成別の最良条件

| 構成 | profile | 条件 | MCC | SPS | CBP F1 | SCI | wall秒 | RSS MiB |
|---|---|---|---:|---:|---:|---:|---:|---:|
| CONTRAfold x CONTRAlign (Alifold) | RNAalifold | `cc_alifold_fth03` | 0.6997 | 0.7707 | 0.5941 | 0.6916 | 1.73 | 13.3 |
| CONTRAfold x CONTRAlign (no Alifold) | off | `cc_noalifold_w2` | 0.6877 | 0.7726 | 0.5844 | 0.6897 | 0.97 | 10.8 |
| ViennaRNA x ProbCons (Alifold) | RNAalifold | `vp_alifold_w2` | 0.6820 | 0.7654 | 0.5778 | 0.6762 | 1.61 | 11.0 |
| ViennaRNA x ProbCons (no Alifold) | off | `vp_noalifold_dd_0p15` | 0.6431 | 0.7633 | 0.5335 | 0.6835 | 0.84 | 10.2 |
| lpc x LinearAlign | LinearAlifold | `cl_linear_dd_0_alifold` | 0.6926 | 0.7787 | 0.5927 | 0.7000 | 5.71 | 20.4 |
| lpc x LinearAlign | off | `cl_linear_ath02` | 0.6881 | 0.7569 | 0.5742 | 0.6721 | 2.47 | 17.4 |
| lpc x LinearCONTRAlign | LinearAlifold | `lpc_contralign_w2_alifold` | 0.6927 | 0.7550 | 0.5833 | 0.6749 | 5.90 | 20.8 |
| lpc x LinearCONTRAlign | off | `lpc_contralign_w8` | 0.6864 | 0.7497 | 0.5789 | 0.6788 | 2.37 | 18.0 |
| lpv x LinearAlign | LinearAlifold | `lpv_linear_dd_0p01_alifold` | 0.6809 | 0.7769 | 0.5801 | 0.6982 | 5.78 | 20.9 |
| lpv x LinearAlign | off | `lpv_linear_fpct05` | 0.6694 | 0.7728 | 0.5690 | 0.6809 | 2.65 | 18.7 |
| lpv x LinearProbCons | LinearAlifold | `lpv_probcons_w2_alifold` | 0.6860 | 0.7536 | 0.5696 | 0.6789 | 5.85 | 20.5 |
| lpv x LinearProbCons | off | `lpv_probcons_dd_0p075` | 0.6648 | 0.7671 | 0.5607 | 0.6760 | 2.32 | 17.5 |

## phase 別の MCC 最良条件

各セルは条件 ID と85ケース平均MCC。

| 構成 | profile | base | parameter | DD RIBOSUM | 最終 RIBOSUM |
|---|---|---|---|---|---|
| CONTRAfold x CONTRAlign (Alifold) | RNAalifold | `cc_alifold_base` (0.6886) | `cc_alifold_fth03` (0.6997) | `cc_alifold_dd_0p05` (0.6927) | `cc_alifold_final_0p005` (0.6885) |
| CONTRAfold x CONTRAlign (no Alifold) | off | `cc_noalifold_base` (0.6811) | `cc_noalifold_w2` (0.6877) | `cc_noalifold_dd_0p075` (0.6826) | `cc_noalifold_final_0p01` (0.6790) |
| ViennaRNA x ProbCons (Alifold) | RNAalifold | `vp_alifold_base` (0.6792) | `vp_alifold_w2` (0.6820) | `vp_alifold_dd_0p02` (0.6792) | `vp_alifold_final_0p005` (0.6791) |
| ViennaRNA x ProbCons (no Alifold) | off | `vp_noalifold_base` (0.6391) | `vp_noalifold_ath005` (0.6425) | `vp_noalifold_dd_0p15` (0.6431) | `vp_noalifold_final_0p005` (0.6373) |
| lpc x LinearAlign | LinearAlifold | `cl_linear_base_alifold` (0.6821) | `cl_linear_w2_alifold` (0.6923) | `cl_linear_dd_0_alifold` (0.6926) | `cl_linear_final_0p005_alifold` (0.6868) |
| lpc x LinearAlign | off | `cl_linear_base` (0.6812) | `cl_linear_ath02` (0.6881) | `cl_linear_dd_0p3` (0.6855) | `cl_linear_final_0p01` (0.6812) |
| lpc x LinearCONTRAlign | LinearAlifold | `lpc_contralign_base_alifold` (0.6803) | `lpc_contralign_w2_alifold` (0.6927) | `lpc_contralign_dd_0p5_alifold` (0.6843) | `lpc_contralign_final_0p005_alifold` (0.6803) |
| lpc x LinearCONTRAlign | off | `lpc_contralign_base` (0.6833) | `lpc_contralign_w8` (0.6864) | `lpc_contralign_dd_0p5` (0.6800) | `lpc_contralign_final_0p01` (0.6846) |
| lpv x LinearAlign | LinearAlifold | `lpv_linear_base_alifold` (0.6611) | `lpv_linear_apct0_alifold` (0.6807) | `lpv_linear_dd_0p01_alifold` (0.6809) | `lpv_linear_final_0p005_alifold` (0.6598) |
| lpv x LinearAlign | off | `lpv_linear_base` (0.6639) | `lpv_linear_fpct05` (0.6694) | `lpv_linear_dd_0p005` (0.6682) | `lpv_linear_final_0p005` (0.6595) |
| lpv x LinearProbCons | LinearAlifold | `lpv_probcons_base_alifold` (0.6679) | `lpv_probcons_w2_alifold` (0.6860) | `lpv_probcons_dd_0p03_alifold` (0.6858) | `lpv_probcons_final_0p005_alifold` (0.6670) |
| lpv x LinearProbCons | off | `lpv_probcons_base` (0.6610) | `lpv_probcons_fpct05` (0.6635) | `lpv_probcons_dd_0p075` (0.6648) | `lpv_probcons_final_0p005` (0.6558) |

## LinearAlifold の対応条件差

linear系の同一条件・同一ケースを対応させた13600組の記述統計。
条件間には共通ケースがあるため、独立標本としての有意性検定は行わない。
非linear系のAlifold有効・無効群は他の引数も異なるため、Alifold単独の効果として扱わない。

| 指標 | 有効−無効の平均差 | 有効が高い組数 |
|---|---:|---:|
| mcc | +0.0110 | 7706/13600 |
| sps | -0.0117 | 4122/13600 |
| sci | +0.0251 | 6978/13600 |
| wall_seconds | +3.4235 | 13401/13600 |
| RSS MiB | +2.7608 | 13595/13600 |

### linear 構成ごとの対応条件差

各構成について条件数×85ケースの対応差を集計した。条件とケースの反復は相関している。

| 構成 | 組数 | ΔMCC | MCC改善組 | ΔSPS | ΔSCI | Δwall秒 | ΔRSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| lpc x LinearAlign | 3400 | +0.0084 | 1969/3400 | -0.0122 | +0.0289 | +3.49 | +2.8 |
| lpc x LinearCONTRAlign | 3400 | +0.0091 | 1954/3400 | -0.0126 | +0.0515 | +3.46 | +2.7 |
| lpv x LinearAlign | 3400 | +0.0081 | 1801/3400 | -0.0134 | +0.0063 | +3.37 | +2.7 |
| lpv x LinearProbCons | 3400 | +0.0182 | 1982/3400 | -0.0084 | +0.0137 | +3.38 | +2.8 |

## 全条件の上位10件

| 順位 | 条件 | MCC | SPS | SCI | wall秒 |
|---:|---|---:|---:|---:|---:|
| 1 | `cc_alifold_fth03` | 0.6997 | 0.7707 | 0.6916 | 1.73 |
| 2 | `lpc_contralign_w2_alifold` | 0.6927 | 0.7550 | 0.6749 | 5.90 |
| 3 | `cc_alifold_dd_0p05` | 0.6927 | 0.7724 | 0.7128 | 1.97 |
| 4 | `cl_linear_dd_0_alifold` | 0.6926 | 0.7787 | 0.7000 | 5.71 |
| 5 | `cc_alifold_dd_0p03` | 0.6926 | 0.7724 | 0.7095 | 1.82 |
| 6 | `cl_linear_dd_0p005_alifold` | 0.6926 | 0.7778 | 0.7025 | 5.74 |
| 7 | `cc_alifold_dd_0p005` | 0.6925 | 0.7736 | 0.7089 | 1.82 |
| 8 | `cc_alifold_dd_0p02` | 0.6925 | 0.7739 | 0.7120 | 1.78 |
| 9 | `cl_linear_dd_0p03_alifold` | 0.6923 | 0.7766 | 0.6969 | 5.75 |
| 10 | `cl_linear_w2_alifold` | 0.6923 | 0.7526 | 0.6716 | 5.97 |

## 最良条件の全引数

### CONTRAfold x CONTRAlign (Alifold) — RNAalifold — `cc_alifold_fth03`

```text
-a CONTRAlign -s CONTRAfold --alifold --dynamic-cbp --ribosum-weight=0.1 --final-ribosum-weight=0 --weight=2 --fold-th=0.3 --max-iter=100
```

### CONTRAfold x CONTRAlign (no Alifold) — off — `cc_noalifold_w2`

```text
-a CONTRAlign -s CONTRAfold --no-alifold --dynamic-cbp --ribosum-weight=0.1 --final-ribosum-weight=0 --fold-pct=0 --weight=2 --max-iter=100
```

### ViennaRNA x ProbCons (Alifold) — RNAalifold — `vp_alifold_w2`

```text
-a ProbCons -s Vienna --alifold --dynamic-cbp --ribosum-weight=0.02 --final-ribosum-weight=0 --fold-th=0.3 --weight=2 --max-iter=100
```

### ViennaRNA x ProbCons (no Alifold) — off — `vp_noalifold_dd_0p15`

```text
-a ProbCons -s Vienna --no-alifold --dynamic-cbp --weight=2 --ribosum-weight=0.15 --final-ribosum-weight=0 --max-iter=100
```

### lpc x LinearAlign — LinearAlifold — `cl_linear_dd_0_alifold`

```text
-a LinearAlign -s lpc --alifold --dynamic-cbp --align-th=0.005 --ribosum-weight=0 --final-ribosum-weight=0 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearAlign — off — `cl_linear_ath02`

```text
-a LinearAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.5 --final-ribosum-weight=0.01 --align-th=0.02 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearCONTRAlign — LinearAlifold — `lpc_contralign_w2_alifold`

```text
-a LinearAlign-CONTRAlign -s lpc --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0.005 --align-pct=0.5 --weight=2 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearCONTRAlign — off — `lpc_contralign_w8`

```text
-a LinearAlign-CONTRAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0.005 --align-pct=0.5 --weight=8 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearAlign — LinearAlifold — `lpv_linear_dd_0p01_alifold`

```text
-a LinearAlign -s lpv --alifold --dynamic-cbp --align-pct=0.5 --ribosum-weight=0.01 --final-ribosum-weight=0 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearAlign — off — `lpv_linear_fpct05`

```text
-a LinearAlign -s lpv --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --align-pct=0.5 --fold-pct=0.5 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearProbCons — LinearAlifold — `lpv_probcons_w2_alifold`

```text
-a LinearAlign-ProbConsRNA -s lpv --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --weight=2 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearProbCons — off — `lpv_probcons_dd_0p075`

```text
-a LinearAlign-ProbConsRNA -s lpv --no-alifold --dynamic-cbp --weight=4 --ribosum-weight=0.075 --final-ribosum-weight=0 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

## 成果物

- 設定: `benchmarks/configs/murlet-dataset1-linearalifold-repair-v2-iter100-beam100.json`
- 実測CSV: `benchmarks/results-murlet-dataset1-linearalifold-repair-v3-iter100-beam100-b5e35428/summary.csv`
