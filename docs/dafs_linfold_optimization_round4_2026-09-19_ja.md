# DAFS LinearFold 最適化 round4（2026-09-19）

## 結論

LinearFold の medium-length 経路に、N-state の三角形 pointer lookup と dense float BPP accumulator を採用した。長さ 64..1024 だけで使い、範囲外は既存の sparse 経路に戻る。N-state の挿入順序と BPP の float 加算順序は維持した。

Murlet dataset1 の85入力×4条件、baseline/candidate 各340 runを fresh に実行した。prediction、quality、非時間 run summary、metrics trace は340/340で一致し、candidate の平均 wall は2.434353秒から2.275618秒へ1.069755倍（6.520636%短縮）だった。正式な結果と340件の比較レコードは [results-linfold-optimization-round4-2026-09-19.json](../benchmarks/results-linfold-optimization-round4-2026-09-19.json) に保存した。

## 固定した実行条件

| 項目 | 値 |
|---|---|
| architecture / compiler | aarch64 / `/usr/bin/c++` GNU 13.3.0 |
| build | Release、`-O3 -DNDEBUG -std=gnu++17` |
| GLPK | 有効、`/home/linuxbrew/.linuxbrew/lib/libglpk.so` |
| ViennaRNA | 2.7.2、`/home/sato-lab.org/satoken/.local/app/viennarna-2.7.2-arm64` |
| OMP | `OMP_NUM_THREADS=1` |
| max-iter / 全 beam | 100 / 100 |
| dataset / 条件 | Murlet dataset1、85入力、LPC/LPV × Alifold有無 |
| 並列度 | baseline/candidateを順番に各8 workers |
| config | original `23a935626952683d9c9bfc51fcd1241b7007a32a9af5189124251410833e581a`、baseline `821f4994e683195692e711112031c82f9d4fa1b4bff0f1fb7b1d4924075d28f0`、candidate `a0cb92788e475e0cb67e1d6aeda3b0bf6c385ee920a620205229b64be96cd303` |

Baseline は `/tmp/dafs-luna-round4-baseline-freeze/dafs`（SHA-256 `ec882997b27178a51a28f778641d6263e29fc2d4ac463bb9e05b51c063dd4ac9`）を凍結して使用した。candidate は `/tmp/dafs-luna-round4-final-live-build/src/dafs`（SHA-256 `c2b62d91354fe58caa7a5c5fd4ea7fa765d61a0553a0a2b98b8a3b4dfcc53b85`）である。fresh raw results は `/tmp/dafs-luna-round4/final-fresh-j8/full-baseline` と `/tmp/dafs-luna-round4/final-fresh-j8/full-candidate` にある。

主要 source hash は次の通りである。

| source | baseline snapshot | final live |
|---|---|---|
| `src/linearfold/fold/linfold.cpp` | `eab97d1e4e2a417f16b88d52472e70adbc495ff1c23889d732fc3b42cd9a0a4a` | `68e4f79845ea787a88f30710727d9c8e8d9a0367938059eeae20dfb5713da921` |
| `src/linearfold/fold/linfold.h` | `661aa479748b72916fd579e1021eba76490122f9b3fffead9a639340add2e5a6` | `011a09ebd41a70266ee7d8c94b47ced5a9ad12c1a07481af92caa8d4ccaf378e` |
| `src/linearfold/param/profile.h` | `36e2c9407ad1db9e45382f9d084f60f5b2c8e0910152a12d9e13fa916dbd7403` | `36e2c9407ad1db9e45382f9d084f60f5b2c8e0910152a12d9e13fa916dbd7403` |
| `src/linfold_wrapper.cpp` | `54303c095eb84752797199dad0e35708f325384140fdc7b65d76f15c7138203b` | `54303c095eb84752797199dad0e35708f325384140fdc7b65d76f15c7138203b` |

`profile.h` は baseline SHA `36e2c940...` に復元した状態を使用した。round4 の production delta は `linfold.{h,cpp}` の bounded lookup と BPP storage である。

## 340 run の比較

数値は85入力の1入力あたり平均。CPU は user+system、RSS は `/usr/bin/time` の `max_rss_kb`、stage は metrics の `stage.seconds` である。`→` は baseline から candidate の順を示す。

| 条件 | wall 秒 | CPU 秒 | RSS KiB | BPP 秒 | progressive 秒 | final consensus 秒 |
|---|---:|---:|---:|---:|---:|---:|
| LPC、Alifoldなし | 1.623412 → 1.540706（1.053681×） | 1.616235 → 1.534706 | 17944.8 → 18087.5 | 0.660491 → 0.582286 | 0.568834 → 0.566032 | 0.000368 → 0.000368 |
| LPC、Alifoldあり | 3.085059 → 2.864000（1.077185×） | 3.078235 → 2.856471 | 21730.0 → 22323.8 | 0.661715 → 0.585044 | 1.923856 → 1.784629 | 0.104403 → 0.099316 |
| LPV、Alifoldなし | 1.812118 → 1.716706（1.055578×） | 1.804824 → 1.710118 | 19338.5 → 19684.4 | 0.681568 → 0.593718 | 0.724634 → 0.721249 | 0.000291 → 0.000291 |
| LPV、Alifoldあり | 3.216824 → 2.981059（1.079088×） | 3.209647 → 2.973647 | 22782.0 → 23341.6 | 0.680353 → 0.594626 | 2.026582 → 1.883521 | 0.104536 → 0.098886 |

全340 runで baseline/candidate ともに status success、metrics validだった。prediction bytes、quality dictionary（MCC/SCI/SPSを含む）、run summary の非時間項目、正規化した metrics trace はすべて340/340一致した。quality の平均差はMCC、SCI、SPSとも全条件で0.0である。時間項目は比較対象から除外した。

## 直接 BPP probe

`/tmp/dafs-luna-round4-final-probe-compare.json` の26 paired rowsを正式な直接 probe とした。profile length 24/63/64/1024/1025、reuse、sequence/gap/constraint変更、LPC、LPVを含む。

- logZ は全行 bitwise 一致。
- BPP pair key set と全 pair float bits は全行一致。
- 最大 absolute / relative BPP difference はともに `0.0`。
- cutoff `0.001, 0.01, 0.05, 0.1` の crossing は0件。
- baseline/candidate JSONL は `/tmp/dafs-luna-round4-final-probe-baseline.jsonl` と `/tmp/dafs-luna-round4-final-probe-candidate.jsonl` に保持した。

## Tests

`ctest --test-dir /tmp/dafs-luna-round4-final-live-build --output-on-failure -j1` は13/13 passした。`LastTest.log` は `/tmp/dafs-luna-round4-final-live-build/Testing/Temporary/LastTest.log`（SHA-256 `e6a065a493c6aa8976e70e30f97e5639f7e0f5ed95cc0c9be71913c864d8f7c5`）に保存した。

## 採用しなかった候補

- 3 energy cache は `/tmp/dafs-luna-round4-rf00017-paired-v2/paired-results.json` と `comparison.json` で確認した。BPP-only に対し LPC は10.107秒→10.120秒、LPVは10.387秒→10.367秒で、改善はノイズ範囲、RSSは約3.5 MiB増だったため棄却した。
- inline callback cache は完全な正式測定がないため、今回の candidateには採用せず性能主張をしていない。
- N-alpha/N-beta を別々に持つ variant は `/tmp/dafs-luna-round4-attribution/comparison.json` と `results.json` で、比較自体は6/6 exactだったが mean RSS +5.91%だった。production は一つの triangular pointer lookupを共有し、重複N-tableを持たない。

## メモリ上の注意

長さ1024では理論上、pair lookup約1 MiB、position-filter suffix約2 MiB、unified N pointer triangular table約4.0 MiB、dense BPP float table約4.0 MiB、profile pair-score dense table約4.0 MiBを使用し得る。64未満または1024超ではこれらのbounded tableを使わない。fresh formal run の RSS change は条件別に+0.80%、+2.73%、+1.79%、+2.46%（全体平均+0.401 MiB、+2.01%）だった。RSSはプロセス全体の測定値なので、個別tableの割当量とは解釈しない。

## 保存と cleanup

最終 JSON、fresh results、final build、baseline freeze/source、direct probe、round4 attribution JSONを保持した。stale partial formal resultsとround4の rejected build/source/probe/object中間物を削除した。round1〜round3の成果物は変更していない。
