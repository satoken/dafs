# Murlet dataset1 LinearAlifold 条件スイープ

85ケース、480条件、40800 run。
全条件で `--max-iter=100`、linear系の5種類のbeamを100に固定した。
順位は85ケースのmacro平均MCC、SCI、SPSの辞書式順序で決めた。
全runで `status=success`、`metrics_valid=True`、`quality_complete=True` を確認した。
全runの開始時刻と実行ファイルのパスをmanifestと照合し、バイナリと設定のSHA256も確認した。
評価バイナリ SHA256（RNAalifold共変異互換修正前）: `a54a94b721a4346a47d5b2194d4c436d059bdeaca960adefc407e36514cab58d`。
設定 SHA256: `a4243213dff3595a81c2f01166c8512cd0eb1430f3c891f39f8ad74d34fb2efa`。実行開始: `2026-09-19T06:03:10.554944+00:00`。
同じケースを条件選択と結果報告に使う探索結果であり、独立したテスト集合の精度ではない。

## 構成別の最良条件

| 構成 | profile | 条件 | MCC | SPS | CBP F1 | SCI | wall秒 | RSS MiB |
|---|---|---|---:|---:|---:|---:|---:|---:|
| CONTRAfold x CONTRAlign (Alifold) | RNAalifold | `cc_alifold_fth03` | 0.6997 | 0.7707 | 0.5941 | 0.6916 | 1.88 | 13.3 |
| CONTRAfold x CONTRAlign (no Alifold) | off | `cc_noalifold_w2` | 0.6877 | 0.7726 | 0.5844 | 0.6897 | 1.03 | 10.7 |
| ViennaRNA x ProbCons (Alifold) | RNAalifold | `vp_alifold_w2` | 0.6820 | 0.7654 | 0.5778 | 0.6762 | 1.71 | 11.1 |
| ViennaRNA x ProbCons (no Alifold) | off | `vp_noalifold_dd_0p15` | 0.6431 | 0.7633 | 0.5335 | 0.6835 | 0.86 | 10.2 |
| lpc x LinearAlign | LinearAlifold | `cl_linear_fth03_alifold` | 0.6349 | 0.7397 | 0.5180 | 0.6816 | 7.11 | 21.3 |
| lpc x LinearAlign | off | `cl_linear_ath02` | 0.6881 | 0.7569 | 0.5742 | 0.6721 | 2.65 | 17.5 |
| lpc x LinearCONTRAlign | LinearAlifold | `lpc_contralign_fth03_alifold` | 0.6415 | 0.7400 | 0.5254 | 0.7033 | 7.06 | 21.4 |
| lpc x LinearCONTRAlign | off | `lpc_contralign_w8` | 0.6864 | 0.7497 | 0.5789 | 0.6788 | 2.53 | 18.1 |
| lpv x LinearAlign | LinearAlifold | `lpv_linear_fth03_alifold` | 0.6355 | 0.7504 | 0.5247 | 0.6897 | 7.16 | 21.7 |
| lpv x LinearAlign | off | `lpv_linear_fpct05` | 0.6694 | 0.7728 | 0.5690 | 0.6809 | 2.79 | 18.7 |
| lpv x LinearProbCons | LinearAlifold | `lpv_probcons_fth03_alifold` | 0.6360 | 0.7416 | 0.5145 | 0.7018 | 7.05 | 21.2 |
| lpv x LinearProbCons | off | `lpv_probcons_dd_0p075` | 0.6648 | 0.7671 | 0.5607 | 0.6760 | 2.48 | 17.5 |

## phase 別の MCC 最良条件

各セルは条件 ID と85ケース平均MCC。

| 構成 | profile | base | parameter | DD RIBOSUM | 最終 RIBOSUM |
|---|---|---|---|---|---|
| CONTRAfold x CONTRAlign (Alifold) | RNAalifold | `cc_alifold_base` (0.6886) | `cc_alifold_fth03` (0.6997) | `cc_alifold_dd_0p05` (0.6927) | `cc_alifold_final_0p005` (0.6885) |
| CONTRAfold x CONTRAlign (no Alifold) | off | `cc_noalifold_base` (0.6811) | `cc_noalifold_w2` (0.6877) | `cc_noalifold_dd_0p075` (0.6826) | `cc_noalifold_final_0p01` (0.6790) |
| ViennaRNA x ProbCons (Alifold) | RNAalifold | `vp_alifold_base` (0.6792) | `vp_alifold_w2` (0.6820) | `vp_alifold_dd_0p02` (0.6792) | `vp_alifold_final_0p005` (0.6791) |
| ViennaRNA x ProbCons (no Alifold) | off | `vp_noalifold_base` (0.6391) | `vp_noalifold_ath005` (0.6425) | `vp_noalifold_dd_0p15` (0.6431) | `vp_noalifold_final_0p005` (0.6373) |
| lpc x LinearAlign | LinearAlifold | `cl_linear_base_alifold` (0.6249) | `cl_linear_fth03_alifold` (0.6349) | `cl_linear_dd_0p075_alifold` (0.6343) | `cl_linear_final_0p01_alifold` (0.6249) |
| lpc x LinearAlign | off | `cl_linear_base` (0.6812) | `cl_linear_ath02` (0.6881) | `cl_linear_dd_0p3` (0.6855) | `cl_linear_final_0p01` (0.6812) |
| lpc x LinearCONTRAlign | LinearAlifold | `lpc_contralign_base_alifold` (0.6248) | `lpc_contralign_fth03_alifold` (0.6415) | `lpc_contralign_dd_0p01_alifold` (0.6264) | `lpc_contralign_final_0p005_alifold` (0.6248) |
| lpc x LinearCONTRAlign | off | `lpc_contralign_base` (0.6833) | `lpc_contralign_w8` (0.6864) | `lpc_contralign_dd_0p5` (0.6800) | `lpc_contralign_final_0p01` (0.6846) |
| lpv x LinearAlign | LinearAlifold | `lpv_linear_base_alifold` (0.6239) | `lpv_linear_fth03_alifold` (0.6355) | `lpv_linear_dd_0p1_alifold` (0.6276) | `lpv_linear_final_0p005_alifold` (0.6234) |
| lpv x LinearAlign | off | `lpv_linear_base` (0.6639) | `lpv_linear_fpct05` (0.6694) | `lpv_linear_dd_0p005` (0.6682) | `lpv_linear_final_0p005` (0.6595) |
| lpv x LinearProbCons | LinearAlifold | `lpv_probcons_base_alifold` (0.6182) | `lpv_probcons_fth03_alifold` (0.6360) | `lpv_probcons_dd_0p15_alifold` (0.6259) | `lpv_probcons_final_0p005_alifold` (0.6173) |
| lpv x LinearProbCons | off | `lpv_probcons_base` (0.6610) | `lpv_probcons_fpct05` (0.6635) | `lpv_probcons_dd_0p075` (0.6648) | `lpv_probcons_final_0p005` (0.6558) |

## LinearAlifold の対応条件差

linear系の同一条件・同一ケースを対応させた13600組の記述統計。
条件間には共通ケースがあるため、独立標本としての有意性検定は行わない。
非linear系のAlifold有効・無効群は他の引数も異なるため、Alifold単独の効果として扱わない。

| 指標 | 有効−無効の平均差 | 有効が高い組数 |
|---|---:|---:|
| mcc | -0.0361 | 5584/13600 |
| sps | -0.0213 | 2921/13600 |
| sci | +0.0291 | 7120/13600 |
| wall_seconds | +4.5921 | 13596/13600 |
| RSS MiB | +3.4608 | 13599/13600 |

### linear 構成ごとの対応条件差

各構成について条件数×85ケースの対応差を集計した。条件とケースの反復は相関している。

| 構成 | 組数 | ΔMCC | MCC改善組 | ΔSPS | ΔSCI | Δwall秒 | ΔRSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| lpc x LinearAlign | 3400 | -0.0435 | 1498/3400 | -0.0224 | +0.0305 | +4.61 | +3.5 |
| lpc x LinearCONTRAlign | 3400 | -0.0412 | 1497/3400 | -0.0220 | +0.0559 | +4.55 | +3.4 |
| lpv x LinearAlign | 3400 | -0.0302 | 1292/3400 | -0.0235 | +0.0115 | +4.64 | +3.4 |
| lpv x LinearProbCons | 3400 | -0.0297 | 1297/3400 | -0.0174 | +0.0183 | +4.57 | +3.6 |

## 全条件の上位10件

| 順位 | 条件 | MCC | SPS | SCI | wall秒 |
|---:|---|---:|---:|---:|---:|
| 1 | `cc_alifold_fth03` | 0.6997 | 0.7707 | 0.6916 | 1.88 |
| 2 | `cc_alifold_dd_0p05` | 0.6927 | 0.7724 | 0.7128 | 1.89 |
| 3 | `cc_alifold_dd_0p03` | 0.6926 | 0.7724 | 0.7095 | 1.80 |
| 4 | `cc_alifold_dd_0p005` | 0.6925 | 0.7736 | 0.7089 | 1.83 |
| 5 | `cc_alifold_dd_0p02` | 0.6925 | 0.7739 | 0.7120 | 1.87 |
| 6 | `cc_alifold_dd_0` | 0.6921 | 0.7742 | 0.7150 | 1.88 |
| 7 | `cc_alifold_ath02` | 0.6917 | 0.7681 | 0.7004 | 1.75 |
| 8 | `cc_alifold_dd_0p01` | 0.6913 | 0.7737 | 0.7123 | 1.91 |
| 9 | `cc_alifold_ath005` | 0.6910 | 0.7683 | 0.7107 | 1.86 |
| 10 | `cc_alifold_base` | 0.6886 | 0.7676 | 0.7108 | 1.84 |

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

### lpc x LinearAlign — LinearAlifold — `cl_linear_fth03_alifold`

```text
-a LinearAlign -s lpc --alifold --dynamic-cbp --ribosum-weight=0.5 --final-ribosum-weight=0.01 --align-th=0.005 --fold-th=0.3 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearAlign — off — `cl_linear_ath02`

```text
-a LinearAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.5 --final-ribosum-weight=0.01 --align-th=0.02 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearCONTRAlign — LinearAlifold — `lpc_contralign_fth03_alifold`

```text
-a LinearAlign-CONTRAlign -s lpc --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0.005 --align-pct=0.5 --fold-th=0.3 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearCONTRAlign — off — `lpc_contralign_w8`

```text
-a LinearAlign-CONTRAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0.005 --align-pct=0.5 --weight=8 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearAlign — LinearAlifold — `lpv_linear_fth03_alifold`

```text
-a LinearAlign -s lpv --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --align-pct=0.5 --fold-th=0.3 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearAlign — off — `lpv_linear_fpct05`

```text
-a LinearAlign -s lpv --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --align-pct=0.5 --fold-pct=0.5 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearProbCons — LinearAlifold — `lpv_probcons_fth03_alifold`

```text
-a LinearAlign-ProbConsRNA -s lpv --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --weight=4 --fold-th=0.3 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearProbCons — off — `lpv_probcons_dd_0p075`

```text
-a LinearAlign-ProbConsRNA -s lpv --no-alifold --dynamic-cbp --weight=4 --ribosum-weight=0.075 --final-ribosum-weight=0 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

## 成果物

- 設定: `/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/configs/murlet-dataset1-linearalifold-repair-v2-iter100-beam100.json`
- 実測CSV: `/home/sato-lab.org/satoken/ghq/github.com/satoken/dafs/benchmarks/results-murlet-dataset1-linearalifold-repair-v2-iter100-beam100/summary.csv`
