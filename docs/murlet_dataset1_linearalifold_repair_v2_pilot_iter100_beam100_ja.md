# Murlet dataset1 LinearAlifold 条件スイープ

3ケース、8条件、24 run。
全条件で `--max-iter=100`、linear系の5種類のbeamを100に固定した。
順位は3ケースのmacro平均MCC、SCI、SPSの辞書式順序で決めた。
全runで `status=success`、`metrics_valid=True`、`quality_complete=True` を確認した。
全runの開始時刻と実行ファイルのパスをmanifestと照合し、バイナリと設定のSHA256も確認した。
評価バイナリ SHA256（RNAalifold共変異互換修正前）: `a54a94b721a4346a47d5b2194d4c436d059bdeaca960adefc407e36514cab58d`。
設定 SHA256: `1f901b937dd60112cb6f3ad1f9719720185a8b83f31af1a2a279fe07842c18b0`。実行開始: `2026-09-19T06:02:18.060189+00:00`。
同じケースを条件選択と結果報告に使う探索結果であり、独立したテスト集合の精度ではない。

## 構成別の最良条件

| 構成 | profile | 条件 | MCC | SPS | CBP F1 | SCI | wall秒 | RSS MiB |
|---|---|---|---:|---:|---:|---:|---:|---:|
| lpc x LinearAlign | LinearAlifold | `cl_linear_ath02_alifold` | 0.6246 | 0.5739 | 0.4096 | 0.7094 | 11.22 | 25.9 |
| lpc x LinearAlign | off | `cl_linear_ath02` | 0.6251 | 0.5867 | 0.4433 | 0.7082 | 2.83 | 21.5 |
| lpc x LinearCONTRAlign | LinearAlifold | `lpc_contralign_ath02_alifold` | 0.6017 | 0.5820 | 0.4192 | 0.7961 | 10.71 | 26.3 |
| lpc x LinearCONTRAlign | off | `lpc_contralign_ath02` | 0.6331 | 0.6003 | 0.4866 | 0.7401 | 3.02 | 21.9 |
| lpv x LinearAlign | LinearAlifold | `lpv_linear_fpct05_alifold` | 0.5919 | 0.6003 | 0.4137 | 0.7011 | 10.88 | 26.7 |
| lpv x LinearAlign | off | `lpv_linear_fpct05` | 0.6350 | 0.6134 | 0.4495 | 0.7147 | 3.25 | 22.8 |
| lpv x LinearProbCons | LinearAlifold | `lpv_probcons_fpct05_alifold` | 0.5796 | 0.5545 | 0.3596 | 0.6998 | 10.25 | 25.9 |
| lpv x LinearProbCons | off | `lpv_probcons_fpct05` | 0.6088 | 0.5845 | 0.4107 | 0.6519 | 3.08 | 21.7 |

## phase 別の MCC 最良条件

各セルは条件 ID と3ケース平均MCC。

| 構成 | profile | base | parameter | DD RIBOSUM | 最終 RIBOSUM |
|---|---|---|---|---|---|
| lpc x LinearAlign | LinearAlifold | — | `cl_linear_ath02_alifold` (0.6246) | — | — |
| lpc x LinearAlign | off | — | `cl_linear_ath02` (0.6251) | — | — |
| lpc x LinearCONTRAlign | LinearAlifold | — | `lpc_contralign_ath02_alifold` (0.6017) | — | — |
| lpc x LinearCONTRAlign | off | — | `lpc_contralign_ath02` (0.6331) | — | — |
| lpv x LinearAlign | LinearAlifold | — | `lpv_linear_fpct05_alifold` (0.5919) | — | — |
| lpv x LinearAlign | off | — | `lpv_linear_fpct05` (0.6350) | — | — |
| lpv x LinearProbCons | LinearAlifold | — | `lpv_probcons_fpct05_alifold` (0.5796) | — | — |
| lpv x LinearProbCons | off | — | `lpv_probcons_fpct05` (0.6088) | — | — |

## LinearAlifold の対応条件差

linear系の同一条件・同一ケースを対応させた12組の記述統計。
条件間には共通ケースがあるため、独立標本としての有意性検定は行わない。
非linear系のAlifold有効・無効群は他の引数も異なるため、Alifold単独の効果として扱わない。

| 指標 | 有効−無効の平均差 | 有効が高い組数 |
|---|---:|---:|
| mcc | -0.0261 | 3/12 |
| sps | -0.0186 | 1/12 |
| sci | +0.0228 | 7/12 |
| wall_seconds | +7.7208 | 12/12 |
| RSS MiB | +4.2028 | 12/12 |

### linear 構成ごとの対応条件差

各構成について条件数×3ケースの対応差を集計した。条件とケースの反復は相関している。

| 構成 | 組数 | ΔMCC | MCC改善組 | ΔSPS | ΔSCI | Δwall秒 | ΔRSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| lpc x LinearAlign | 3 | -0.0005 | 2/3 | -0.0128 | +0.0012 | +8.39 | +4.3 |
| lpc x LinearCONTRAlign | 3 | -0.0314 | 0/3 | -0.0183 | +0.0559 | +7.69 | +4.4 |
| lpv x LinearAlign | 3 | -0.0431 | 0/3 | -0.0131 | -0.0137 | +7.63 | +3.9 |
| lpv x LinearProbCons | 3 | -0.0292 | 1/3 | -0.0300 | +0.0479 | +7.18 | +4.2 |

## 全条件の上位10件

| 順位 | 条件 | MCC | SPS | SCI | wall秒 |
|---:|---|---:|---:|---:|---:|
| 1 | `lpv_linear_fpct05` | 0.6350 | 0.6134 | 0.7147 | 3.25 |
| 2 | `lpc_contralign_ath02` | 0.6331 | 0.6003 | 0.7401 | 3.02 |
| 3 | `cl_linear_ath02` | 0.6251 | 0.5867 | 0.7082 | 2.83 |
| 4 | `cl_linear_ath02_alifold` | 0.6246 | 0.5739 | 0.7094 | 11.22 |
| 5 | `lpv_probcons_fpct05` | 0.6088 | 0.5845 | 0.6519 | 3.08 |
| 6 | `lpc_contralign_ath02_alifold` | 0.6017 | 0.5820 | 0.7961 | 10.71 |
| 7 | `lpv_linear_fpct05_alifold` | 0.5919 | 0.6003 | 0.7011 | 10.88 |
| 8 | `lpv_probcons_fpct05_alifold` | 0.5796 | 0.5545 | 0.6998 | 10.25 |

## 最良条件の全引数

### lpc x LinearAlign — LinearAlifold — `cl_linear_ath02_alifold`

```text
-a LinearAlign -s lpc --alifold --dynamic-cbp --ribosum-weight=0.5 --final-ribosum-weight=0.01 --align-th=0.02 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearAlign — off — `cl_linear_ath02`

```text
-a LinearAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.5 --final-ribosum-weight=0.01 --align-th=0.02 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearCONTRAlign — LinearAlifold — `lpc_contralign_ath02_alifold`

```text
-a LinearAlign-CONTRAlign -s lpc --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0.005 --align-pct=0.5 --align-th=0.02 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpc x LinearCONTRAlign — off — `lpc_contralign_ath02`

```text
-a LinearAlign-CONTRAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0.005 --align-pct=0.5 --align-th=0.02 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearAlign — LinearAlifold — `lpv_linear_fpct05_alifold`

```text
-a LinearAlign -s lpv --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --align-pct=0.5 --fold-pct=0.5 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearAlign — off — `lpv_linear_fpct05`

```text
-a LinearAlign -s lpv --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --align-pct=0.5 --fold-pct=0.5 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearProbCons — LinearAlifold — `lpv_probcons_fpct05_alifold`

```text
-a LinearAlign-ProbConsRNA -s lpv --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --weight=4 --fold-pct=0.5 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

### lpv x LinearProbCons — off — `lpv_probcons_fpct05`

```text
-a LinearAlign-ProbConsRNA -s lpv --no-alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --weight=4 --fold-pct=0.5 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100
```

## 成果物

- 設定: `benchmarks/configs/murlet-dataset1-linearalifold-repair-v2-pilot-iter100-beam100.json`
- 実測CSV: `benchmarks/results-murlet-dataset1-linearalifold-repair-v2-pilot-iter100-beam100/summary.csv`
