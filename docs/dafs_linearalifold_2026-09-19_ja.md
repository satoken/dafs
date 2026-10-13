# DAFS の linear profile folding 修正

## 実装

`-s lpv` / `-s lpc` と `--alifold` の併用により、最終構造予測に使う
プロフィールの分配関数と塩基対確率を LinFold の inside/outside で計算する。
`--linfold-beam` がプロフィール探索にも適用される。`--no-alifold` は既定値。
非 linear 系では従来の RNAalifold を使う。

プロフィールは両モデルで共通の Turner 熱力学モデルを使う。各配列へ
ループを射影し、エネルギーの配列平均を `-E/RT` に変換する（37℃）。
ViennaRNA 2.7 以降の最近接エネルギー評価関数を利用するが、プロフィールの
探索、枝刈り、分配関数、outside、周辺確率は LinFold で計算する。
ViennaRNA の RNAalifold 動的計画法を呼び出す実装ではない。

列ペアの候補には canonical な塩基対を持つ配列が最低1本必要であり、
`2 × incompatible + double-gap < 配列数` を要求する。非 canonical な
行のエネルギー評価は Turner の type 7 を使う。内部ループの探索範囲は
列数ではなく、ギャップを除いた長さの配列平均で制限する。
投影後に3塩基未満となるヘアピンには6 kcal/molを与える。
3塩基ヘアピンには一般的な triloop の規則を適用し、配列特異的な
tri/tetra/hexa-loop ボーナスは使わない。

塩基対ごとに ViennaRNA RNAalifold と同じ既定の共変異項と、非対応塩基・double-gap
へのペナルティを加える。canonical pair 型間の Hamming 距離を配列数で正規化し、
RNAalifold の PF と同じ centikcal スケールへ変換する。DAFS の RIBOSUM85-60 は
alignment と最終復号の目的関数にだけ使い、プロフィール分配関数には加えない。
公開 LinearAlifold の適応的 RIBOSUM / BL* と同一の数値モデルではない。
LazyOutside は使わず、枝刈り後の探索グラフを逆向きにたどる。

## 修正した不具合

- outside/BPP の map 読み出しで、inside 時に枝刈りされた状態を再生成していた。
  挿入しない参照に変更し、繰り返し outside を実行する前に beta を初期化する。
- 長いヘリックスのループ変数が shadowing により未初期化になり、
  左座標の減算には unsigned のアンダーフローもあった。
- 確率の事後的な総和正規化を除去し、分配関数から得られる周辺確率を返す。
  非有限値や確率質量の過大値は例外として検出する。
- 外側の先頭1〜2塩基の初期化に制約を適用する。
- 固定塩基対の左右入れ替えで参照が alias する問題を修正する。
- DAFS の構造制約生成で、右括弧を再処理して固定塩基対を消していた。
- `--alifold --bp-update1` で、混合確率から予測された少数配列支持の塩基対を
  プロフィールに強制すると失敗していた。DAFS 内部の反復更新に限り、
  モデルが許容しない予測塩基対を非対合として扱う。各配列の更新は保持する。
  wrapper を直接呼ぶ際の固定塩基対検証は既定で厳密なまま。
- ギャップの累積座標、隣接塩基、非 canonical 行のエネルギー、RT による
  無次元化を明示的に実装する。

`calculate_profile` は logZ を返す。プロフィールの確率出力しきい値は
wrapper の設定に従う。構造制約は `?`、`.`、`(`、`)` を受け付け、
不正・実現不能な固定塩基対を黙って解除しない。
上記の DAFS 内部更新は `relax_unsupported_pairs=true` を明示する。

## 検証

- 単一行プロフィール130配列について、独立した ViennaRNA の分配関数
  DP と logZ・全塩基対確率を照合する。同じ熱力学モデルを使うため、
  検証側では特殊ループのモチーフを空にし、PF smoothing を無効化する。
- 2構造だけの解析的なアンサンブル、全非対合制約、固定塩基対、
  34塩基対ヘリックス、全ギャップ列、45列にわたるギャップを含む内部ループ、
  共変異プロフィールを確認する。
- 2行の解析的アンサンブルで、共変異の正規化、非 canonical 行、両端ギャップ行
  の logZ と BPP を照合する。少数配列支持の塩基対の厳密／内部更新用制約も確認する。
- linear `--alifold` と `--bp-update1` のパイプラインを確認する。
- メモリ・未定義動作の検出には AddressSanitizer / UndefinedBehaviorSanitizer を使う。
- Luna Max が作成した別の検証では、長さ0〜12の非交差構造を全列挙し、
  `max_helix=0,1,2,3,30` の分配関数と全塩基対確率を照合する。
  非対合制約・列ペアのフィルタも確認し、beam `0,1,2,3,5,100` では
  logZ の有限差分と周辺確率から計算した期待値を比較する。
  これらを含む CTest 13件が成功した。上記の beam の変更は単体検証のみで、
  精度・速度評価の beam は100に固定する。

### RNAalifold との短配列比較

ViennaRNA 2.7.2 の `vrna_pf_alifold`（既定モデル）と、LPV の
`calculate_profile`（beam=100）を同じ短いアラインメントで直接比較した。
`logZ` は `-G/RT`、BPP は全列ペアの絶対差で集計した。単一行ではほぼ一致するが、
複数行でも保存 pair と compensatory pair はほぼ一致する。

| 配列プロフィール | ΔlogZ (Linear−RNAalifold) | BPP RMS差 | 構造（BPP>0.5） |
|---|---:|---:|---|
| `GGGAAACCC`（1行） | −0.00056 | 0.00016 | `(((...)))` |
| `GGGAAACCC` × 2 | −0.00056 | 0.00016 | 一致 |
| `GGGAAACCC` / `GGGAAUCCC` | −0.00070 | 0.00021 | 一致 |
| `GGGAAACCC` / `CCCAAAGGG` | −0.00024 | 0.00008 | 一致 |
| `GGGAAACCC` / `GCGAAAGGC` / `GUGAAACAC` | +0.21315 | 0.01950 | RNAalifold `.........`、Linear `((.....))` |

この比較で確認した差の大部分は、変更前にプロフィールへ RIBOSUM 自己スコアを
加えていたことが原因だった。現在は RNAalifold の既定共変異式へ合わせ、保存 pair
と compensatory pair の短配列ケースが同じ logZ/BPP になることを確認した。
ViennaRNA の特殊ループ処理、ギャップを含む投影、枝刈り後の計算順序はなお LinearFold
側の実装であり、任意の長いプロフィールでのビット単位互換を主張するものではない。

## 評価設計

先行タスク「murlet実験のMCC低下原因を調査」に対応するローカル資料
`docs/murlet_dataset1_condition_sweep_2026-09-18_ja.md` と320条件設定を参照した。
タスクのメッセージ本文は取得できなかったため、評価条件の根拠には
リポジトリ内の設定・結果資料を使った。

Murlet dataset1 の85ケース、既存320条件を `--max-iter=100` に統一し、
linear 160条件に同じ引数の `--alifold` 有効版を追加する。
合計480条件、40,800 run。linear 系の `--align-beam`、`--align-dd-beam`、
`--linfold-beam`、`--fold-dd-beam`、`--fold-final-beam` はすべて100。
各条件・ケースは1回、`OMP_NUM_THREADS=1`、SCI採点は ViennaRNA 2.7.2。

評価担当は GPT-5.6 Luna / max。修正前の途中結果は採用せず、別の出力先で
24 run の予備評価から再実行する。実行ファイル・設定・入力・参照・採点ツール
のハッシュを記録し、再開時に異なる実行ファイルの結果を再利用しない。
全件の成功とケース被覆を検証した後、85ケース macro 平均 MCC、SCI、SPS
の順で順位を決める。条件選択と報告に同じケースを使う探索評価であり、
独立したテスト集合の精度ではない。

互換修正前の全件評価バイナリ SHA256 は
`a54a94b721a4346a47d5b2194d4c436d059bdeaca960adefc407e36514cab58d`。
今回の RNAalifold 共変異互換修正後の検証ビルド SHA256 は
`b5e3542853b6e368a88d483787365745d1805adedf36ebcb2a74e8ddbbd6fee7`。
C/C++ ソースとビルド設定109ファイルのハッシュは
`benchmarks/linearalifold-repair-v2-source-sha256.json` に保存した。
互換修正前バイナリで予備評価24件が成功した。これは条件設計と runner の
動作確認用であり、今回の共変異互換修正後の精度値としては扱わない。今回の修正後は
短配列の直接比較と CTest 13/13 を別に実施した。設定・既存結果は以下に保存している。

- 全件設定: `benchmarks/configs/murlet-dataset1-linearalifold-repair-v2-iter100-beam100.json`
- 予備評価: [24 run の結果](murlet_dataset1_linearalifold_repair_v2_pilot_iter100_beam100_ja.md)
- 全件出力: `benchmarks/results-murlet-dataset1-linearalifold-repair-v2-iter100-beam100/`

## 互換修正前の全件評価結果

互換修正前の40,800 runはすべて成功した。全480条件に85ケースがそろい、
`metrics_valid` と `quality_complete` は全件 true。
CSVの全列を各 `result.json` と照合し、実際の入力・参照ファイル、
実行バイナリ、設定、SCI採点ツールのハッシュを検証した。
評価終了後も、バイナリおよび109個のソース・ビルド関連ファイルに変更はなかった。
この結果は今回の RNAalifold 共変異互換修正前の記録であり、今回の修正後の
精度評価には再利用しない。修正後の CTest は13/13成功した。

| 区分 | 最良条件 | 85ケース平均MCC |
|---|---|---:|
| 全条件 | `cc_alifold_fth03` | 0.6997 |
| linear・LinearAlifold有効 | `lpc_contralign_fth03_alifold` | 0.6415 |
| linear・無効 | `cl_linear_ath02` | 0.6881 |

上の行間では選択された引数が異なる。alifold 以外をそろえた13,600組の
対応比較では、有効−無効の平均差は MCC **−0.0361**、SPS **−0.0213**、
SCI **＋0.0291**、wall time **＋4.5921秒/ケース**、RSS **＋3.4608 MiB** だった。
MCC が改善したのは5,584/13,600組。今回の条件範囲では、4つのlinear構成すべてで
平均MCCが低下した。LinearAlifold は利用可能になったが、今回の結果から
一律の精度向上は支持されない。既定の `--no-alifold` は変更していない。

構成別・phase別の結果、対応差、最良条件の全引数は
[全件評価レポート](murlet_dataset1_linearalifold_repair_v2_iter100_beam100_ja.md) に記録した。
最終集計は `python3 benchmarks/summarize_linearalifold_sweep.py` で再生成できる。

評価したビルドは `/tmp/dafs-linearalifold-repair-build`。
ViennaRNA 2.7.2 の pkg-config 設定を指定し、CMake の Release / BUILD_TESTING=ON
でビルドした。評価の実行コマンドは以下。

```sh
python3 benchmarks/run.py \
  --config benchmarks/configs/murlet-dataset1-linearalifold-repair-v2-iter100-beam100.json \
  --jobs 20
```

この設定に記録されたバイナリは上記の評価ビルドを指す。別の場所で再ビルドした
場合は `--binary <dafsのパス>` を明示する。実行ファイルのハッシュが変われば、
既存の計算結果は再利用されない。

## 参考

- [LinearAlifold 論文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11377157/)
- [LinearAlifold 公開実装](https://github.com/LinearFold/LinearAlifold)
  （照合 commit: `138feef3e04e18d50a306ed79d4e9dc2ba28f235`）
