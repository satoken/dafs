# LinearPartition / LinearAlifold 高速化 round3（2026-09-19）

## 結論

Profile のエネルギー計算の簡略化と LinearAlifold の profile pair-score dense cache を採用する。正式な比較は GLPK 有効の同一構成で実施し、85入力×4条件の340 runを比較した。

予測FASTA、品質値、DD trace と run summary の非時間項目、直接測定した profile BPP はすべて baseline と一致した。Alifold 条件の wall は LPC で 1.066946倍、LPV で 1.063399倍だった。

機械可読な全結果は [results-linfold-optimization-round3-2026-09-19.json](../benchmarks/results-linfold-optimization-round3-2026-09-19.json) に保存した。

## 対象と変更

今回の対象は LinearPartition と LinearAlifold の確率計算（PF）である。LinFold core の `linfold.cpp` は開始時の内容と同じである。

- `profile.h` の profile 内部ループ、pair type、平均化の fastpath を使う。平均化の除算と signed zero は元の計算結果を保つ形にした。
- `linfold_wrapper.cpp` の LinearAlifold profile pair-score を、長さ64～1024では dense cache に保持する。範囲外は既存の map 経路を使う。
- profile fastpath の直接 probe は、logZ と全 BPP の bitwise 一致を確認した。

## 固定した build と入力

baseline は現行 source から先に binary を凍結した。正式 candidate は同じ live tree を次の構成で再 build した。

| 項目 | 値 |
|---|---|
| architecture | ARM64 (`aarch64`) |
| compiler | `/usr/bin/c++` GNU 13.3.0 |
| CMake | Release, `-O3 -DNDEBUG -std=gnu++17` |
| ViennaRNA | 2.7.2 |
| GLPK | 有効（`/home/linuxbrew/.linuxbrew/lib/libglpk.so`） |
| `OMP_NUM_THREADS` | 1 |
| max-iter | 100 |
| 全 beam | 100 |
| dataset | Murlet dataset1、85入力 |
| 条件 | LPC/LPV × Alifoldあり/なし |

設定ファイルの SHA-256 は `23a935626952683d9c9bfc51fcd1241b7007a32a9af5189124251410833e581a` である。formal run は8並列で、各入力を1回測った。

### Binary と source hash

| 対象 | path | SHA-256 |
|---|---|---|
| baseline | `/tmp/dafs-luna-round3-baseline.LW7gSD/dafs` | `cd3e20c9d7267440bca9440e2d6b61748117dcdb0552e1398b180ff662654f19` |
| formal candidate | `/tmp/dafs-linfold-fast-build/src/dafs` | `ec882997b27178a51a28f778641d6263e29fc2d4ac463bb9e05b51c063dd4ac9` |
| `src/linearfold/fold/linfold.cpp` | live/final source | `eab97d1e4e2a417f16b88d52472e70adbc495ff1c23889d732fc3b42cd9a0a4a` |
| `src/linfold_wrapper.cpp` | live/final source | `54303c095eb84752797199dad0e35708f325384140fdc7b65d76f15c7138203b` |
| `src/linearfold/param/profile.h` | live/final source | `36e2c9407ad1db9e45382f9d084f60f5b2c8e0910152a12d9e13fa916dbd7403` |

## Formal 85×4 結果

時間は各条件の1入力あたり平均秒である。`wall` と `CPU` は DAFS process、`DD` は metrics の DD 時間、`BPP` は初期 base-pair probability stage、`progressive` と `final fold` は metrics の stage 時間である。

| 条件 | wall 秒（baseline→candidate） | wall speedup | CPU 秒（baseline→candidate） | DD 秒（baseline→candidate） | BPP 秒（baseline→candidate） | progressive 秒（baseline→candidate） | final fold 秒（baseline→candidate） |
|---|---:|---:|---:|---:|---:|---:|---:|
| LPC、Alifoldなし | 1.648118→1.643647 | 1.002720× | 1.641294→1.637294 | 0.572664→0.570944 | 0.674520→0.670430 | 0.575626→0.573817 | 0.000371→0.000376 |
| LPC、Alifoldあり | 3.320588→3.112235 | 1.066946× | 3.313765→3.106118 | 0.539572→0.526687 | 0.668864→0.666734 | 2.132990→1.941216 | 0.120069→0.105861 |
| LPV、Alifoldなし | 1.823176→1.825529 | 0.998711× | 1.816471→1.818588 | 0.726544→0.723371 | 0.684978→0.688654 | 0.731705→0.728644 | 0.000293→0.000297 |
| LPV、Alifoldあり | 3.453294→3.247412 | 1.063399× | 3.446706→3.240353 | 0.642977→0.624663 | 0.686368→0.690180 | 2.242313→2.043568 | 0.118641→0.105134 |

LPC/LPV の Alifoldあり条件では progressive alignment と final consensus folding の短縮が wall 短縮に寄与した。Alifoldなしの LPV は単発測定の揺らぎ範囲で、wall は0.13%増だった。

### 予測・品質・trace

formal comparison は次の件数だった。

| 項目 | 一致件数 |
|---|---:|
| 実行、paired | 340/340 |
| prediction FASTA bytes | 340/340 |
| quality dict（MCC/SCI/SPSを含む） | 340/340 |
| run summary の非時間項目 | 340/340 |
| DD trace events | 340/340 |
| metrics errors / missing / unexpected | 0 / 0 / 0 |

per-case の MCC、SCI、SPS 最大絶対差はすべて `0.0` だった。品質平均の差も全4条件で `0.0` だった。

### 直接 BPP probe

RF00017 の aligned profile を `ProfileNearestNeighbor` に入力し、baseline と formal candidate の BPP を probability cutoff 0 で比較した。

- pair entries: 7537 / 7537、pair key set 一致
- logZ: `132.42389462173659` / `132.42389462173659`、bitwise 一致
- 最大絶対 BPP 差: `0.0`
- 最大相対 BPP 差: `0.0`
- cutoff 0.001、0.01、0.05、0.1 の跨ぎ: すべて0件

## Tests

GLPK 有効の `/tmp/dafs-linfold-fast-build` で CTest 13/13 を通過した。`linear_models_fail_fast`、`linear_profile_oracle`、`linear_partition_exhaustive`、`relaxed_bounds`、`gradient_manager`、`linear_nussinov`、`ribosum_pair_match`、benchmark scorer/runner の3件、`reject_negative_max_iter`、`linear_pipeline_smoke`、`zero_max_iter_uses_ip` を含む。

## 候補の判定

- A（全 model の node posterior 集約）は棄却した。85×4で prediction/quality は325/340、最悪は `lpv_linear_fpct05_alifold/murlet-dataset1-RF00230-2/rep-01` の MCC `0.3332227733→0.2139047960`（差 `-0.1193179773`）だった。
- A（Profile 限定）も pilot の `cl_linear_ath02_alifold/murlet-dataset1-RF00017-0/rep-01` で MCC `0.5075759900→0.4781138791`（差 `-0.0294621109`）となり棄却した。
- C（outside/BPP で beta=-inf になる到達不能先のエネルギー計算を省略）は pilot 24/24で完全一致したが、wall の改善が明確でないため今回の candidate から外した。
- E（`allowed_pairs` bitset の連続配置）は Options のサイズ増加（208→8376 bytes）に対する production 効果が小さいため保留した。
- B（inside/outside/BPP の制約・table 準備を共有）は gprof の self-sample 根拠が小さく、今回実装していない。

A/C の pilot と rejected A の full は探索用 build（GLPKなし）で実行した。正式な速度と採用判定は、GLPK を baseline と揃えた formal candidate の結果だけに基づく。

## Profiling

`perf` は環境の `perf_event_paranoid=4` により使用できなかった。gprof は A-like 探索候補の hotspot の形を確認するために使い、baseline や formal candidate の速度値には使っていない。サンプル比率は profile callback 13.65%、Profile outside 9.75%、Profile add-alpha 9.01%、decode 8.36%、Profile inside 5.48%、logsumexp 3.71% だった。

raw run、比較 CSV/JSON、manifest は `/tmp/dafs-luna-round3-baseline-85`、`/tmp/dafs-luna-round3-profileD-glpk-85`、`/tmp/dafs-luna-round3-profileD-glpk-85-compare` に保存した。

結果保存後、今回作成した隔離 build/source、棄却候補の raw run、旧 normalization の中断 run、temporary BPP probe binary/source を削除した。正式 raw、pilot raw、直接 BPP 出力、baseline frozen binary、formal build は保持している。
