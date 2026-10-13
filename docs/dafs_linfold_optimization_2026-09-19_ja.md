# LinearPartition / LinearAlifold 高速化（2026-09-19）

## 変更

- 両モジュールが共有する `LinFold` の beam 枝刈りを全ソートから `nth_element` に変更し、候補バッファを再利用した。スコア降順・位置昇順の同点優先順を維持する。
- `LinFold` の inside、outside、BPP、Viterbi で繰り返す塩基対可否判定を、配列長 64～1024 の場合に一度計算して保持する。最大 `(1024+1)^2` byte の表を使い、範囲外は従来の判定へ戻る。
- LinearAlifold の複数行プロファイルでも、長さ1024以下の塩基対可否を保持する。プロファイルの行や制約が変わるたびに表を作り直す。
- `LinFoldWrapper` は折りたたみエンジンと枝刈り用バッファを再利用する。新しい配列・プロファイルではパラメータモデルを更新し、以前の DP 状態を消去する。

## 検証方法

比較対象は変更前の frozen binary（SHA-256 `b5e3542853b6e368a88d483787365745d1805adedf36ebcb2a74e8ddbbd6fee7`）と変更後の binary（SHA-256 `e7444644aae21b973a98d261d7defc47f73dbfa0193feb361a124883e5554070`）。ARM64、ViennaRNA 2.7.2、`OMP_NUM_THREADS=1`、8並列で実行した。

Murlet dataset1 の85件について、LPC・LPVそれぞれ Alifold あり／なしの4条件を1回ずつ測った。各 binary 340 run、合計680 run。全条件で `max-iter=100`、`align-beam=100`、`align-dd-beam=100`、`linfold-beam=100`、`fold-dd-beam=100`、`fold-final-beam=100` を固定した。[実行設定](../benchmarks/configs/murlet-dataset1-linfold-optimization-iter100-beam100.json)のSHA-256は `23a935626952683d9c9bfc51fcd1241b7007a32a9af5189124251410833e581a`。

予測FASTA、品質値、時刻項目を除くDD trace、および時刻・実行ファイル識別子を除くrun summaryは、変更前後で **340/340 件すべて一致**。全runで成功、metrics valid、品質評価完了。加えて `ctest` は13/13通過し、長さ64・80の単一行プロファイルについて分配関数とBPPをViennaRNAの非beam計算と比較した。

## 85件平均

`wall` と `CPU` はDAFSプロセスの秒数で、外部SCI評価時間を含まない。速度比は変更前÷変更後。初期BPP段階がLinearPartition、最終構造段階のAlifoldありはLinearAlifoldの処理を含む。

| 条件 | wall 秒（前→後） | wall 速度比 | CPU 秒（前→後） | 初期BPP 秒（前→後） | 最終構造 秒（前→後） | 最大RSS MiB（前→後） |
|---|---:|---:|---:|---:|---:|---:|
| LPC、Alifoldなし | 1.735→1.638 | 1.059× | 1.730→1.632 | 0.773→0.675（1.146×） | — | 17.52→17.55 |
| LPC、Alifoldあり | 4.886→3.580 | 1.365× | 4.879→3.573 | 0.773→0.673（1.148×） | 0.700→0.167（4.186×） | 20.39→21.12 |
| LPV、Alifoldなし | 1.909→1.820 | 1.048× | 1.902→1.813 | 0.785→0.695（1.129×） | — | 18.74→18.91 |
| LPV、Alifoldあり | 4.922→3.718 | 1.324× | 4.914→3.711 | 0.783→0.697（1.125×） | 0.626→0.162（3.876×） | 21.09→22.06 |

Alifoldありの漸進アラインメント段階もLPCで3.019→2.348秒（1.286×）、LPVで3.109→2.455秒（1.266×）。プロファイル分配関数の再計算がこの段階に含まれる。

条件ごとの全測定値と一致判定は[機械可読結果](../benchmarks/results-linfold-optimization-2026-09-19.json)に保存した。8並列での1回測定であり、速度比はこの入力群・設定・マシンに対する観測値である。

再実行には次のコマンドを使う。`--output-dir` は既存結果と重ならない場所を指定する。

```sh
python3 benchmarks/run.py --config benchmarks/configs/murlet-dataset1-linfold-optimization-iter100-beam100.json --binary /tmp/dafs-linearalifold-repair-build/src/dafs --output-dir /tmp/dafs-linfold-rerun-baseline --jobs 8
python3 benchmarks/run.py --config benchmarks/configs/murlet-dataset1-linfold-optimization-iter100-beam100.json --binary /tmp/dafs-linfold-fast-build/src/dafs --output-dir /tmp/dafs-linfold-rerun-final --jobs 8
```
