# N/L × Alifold：Nussinov最適化の速度比較

## 結果

時間はDAFSの `run_summary.seconds` を3回計測した中央値です。タイマーはFASTA読込後に始まり、確率計算・アラインメント・双対分解・最終コンセンサス折りたたみを含みます。速度比は「旧版の時間 ÷ 現行版の時間」で、1を超えると現行版が高速です。

| 長さ / 設定 | 旧版 (秒) | 現行DD 1スレッド (秒) | 速度比 | 現行DD 3スレッド (秒) | 速度比 | DD時間：旧版 → 1スレッド → 3スレッド (秒) | DD速度比 (3スレッド) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 / N0 | 0.070 | 0.078 | 0.90× | 0.069 | 1.02× | 0.025 → 0.033 → 0.024 | 1.07× |
| 100 / N1 | 0.128 | 0.135 | 0.95× | 0.129 | 1.00× | 0.019 → 0.026 → 0.019 | 1.02× |
| 100 / L0 | 0.201 | 0.195 | 1.03× | 0.185 | 1.09× | 0.084 → 0.078 → 0.068 | 1.24× |
| 100 / L1 | 0.395 | 0.391 | 1.01× | 0.386 | 1.02× | 0.063 → 0.059 → 0.054 | 1.18× |
| 250 / N0 | 0.515 | 0.534 | 0.96× | 0.485 | 1.06× | 0.106 → 0.124 → 0.075 | 1.42× |
| 250 / N1 | 1.091 | 1.116 | 0.98× | 1.072 | 1.02× | 0.095 → 0.120 → 0.073 | 1.30× |
| 250 / L0 | 1.327 | 1.266 | 1.05× | 1.168 | 1.14× | 0.407 → 0.350 → 0.252 | 1.62× |
| 250 / L1 | 3.373 | 3.339 | 1.01× | 3.301 | 1.02× | 0.294 → 0.269 → 0.218 | 1.35× |
| 500 / N0 | 2.549 | 2.532 | 1.01× | 2.358 | 1.08× | 0.411 → 0.397 → 0.220 | 1.87× |
| 500 / N1 | 5.322 | 5.326 | 1.00× | 5.170 | 1.03× | 0.377 → 0.387 → 0.217 | 1.74× |
| 500 / L0 | 4.185 | 3.989 | 1.05× | 3.689 | 1.13× | 1.328 → 1.130 → 0.830 | 1.60× |
| 500 / L1 | 11.645 | 11.486 | 1.01× | 11.276 | 1.03× | 0.917 → 0.753 → 0.535 | 1.71× |
| 1000 / N0 | 14.076 | 13.663 | 1.03× | 12.997 | 1.08× | 1.876 → 1.452 → 0.792 | 2.37× |
| 1000 / N1 | 28.582 | 28.222 | 1.01× | 27.618 | 1.03× | 1.844 → 1.429 → 0.792 | 2.33× |
| 1000 / L0 | 10.564 | 10.064 | 1.05× | 9.174 | 1.15× | 3.276 → 2.780 → 1.896 | 1.73× |
| 1000 / L1 | 32.557 | 32.116 | 1.01× | 31.535 | 1.03× | 2.329 → 1.893 → 1.280 | 1.82× |

## 設定

すべて `--allow-isolated --fold-decoder=Nussinov --seed=42 --max-iter=100` を共通指定しました。Nは非線形の `CONTRAlign`、Lは `LinearAlign` を使います。0はAlifoldなし、1はAlifoldありです。使用した推奨パラメータは次の通りです。

- **N0** (`cc_noalifold_w2`): `-a CONTRAlign -s CONTRAfold --no-alifold --dynamic-cbp --ribosum-weight=0.1 --final-ribosum-weight=0 --fold-pct=0 --weight=2`
- **N1** (`cc_alifold_fth03`): `-a CONTRAlign -s CONTRAfold --alifold --dynamic-cbp --ribosum-weight=0.1 --final-ribosum-weight=0 --weight=2 --fold-th=0.3`
- **L0** (`cl_linear_ath02`): `-a LinearAlign -s lpc --no-alifold --dynamic-cbp --ribosum-weight=0.5 --final-ribosum-weight=0.01 --align-th=0.02 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100`
- **L1** (`cl_linear_dd_0_alifold`): `-a LinearAlign -s lpc --alifold --dynamic-cbp --align-th=0.005 --ribosum-weight=0 --final-ribosum-weight=0 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100`

## 計測方法

- `benchmarks/data/rfam15.1-long/inputs/rfam15.1-long-pk-RF00177-n1629.fa` 内のRFAM配列から4配列を選び、先頭100、250、500、1,000 ntを入力にしました。配列長を変えた比較用であり、独立した生物学的評価セットではありません。入力FASTAはローカル生成物としてGit管理対象外です。元データセットの出所とSHA-256は `benchmarks/data/rfam15.1-long/manifest.json` にあります。
- 各長さ・設定を3回計測し、実行順は反復ごとに入れ替えました。全実行をCPU 15–17に固定し、`OMP_NUM_THREADS=3`、`OMP_DYNAMIC=FALSE`、`OMP_PROC_BIND=TRUE`、`OMP_PLACES=cores` を指定しました。
- `before` はNussinovのsparsificationとDD分割へのスレッド割り当てを入れる前の保存バイナリです。`after-1` は現行DAFSに `--dd-threads=1 --fold-dp-threads=1`、`after-3` は `--dd-threads=3 --fold-dp-threads=1` を指定しました。
- 現行版の1/3スレッド間で予測結果と時間以外のDDメトリクスが一致することを実行スクリプトが確認しています。48通りすべて各3回完了し、エラーはありません。
- DAFS内部時間はプロセス起動とFASTA読込を除きます。`DD時間` は progressive alignment の3回の双対分解の合計です。

## バイナリ識別情報

- Before binary (SHA-256 `6511ff1a6ce54ab124248eab5910ec3b13293e632799cba765dadcd6d807a1cc`)
- Current binary built from this working tree (SHA-256 `40ce723e467a17d11da206a4019df281b63457ac19e0543219b2bdb885f557d3`)
