# IPknot Nussinovの同梱とデフォルト化（2026-10-13）

> 同日追記：非線形Nussinovの既定は、旧DAFSの`DAFSNussinov`（`d7f4b07`由来）に戻した。線形`LinearNussinov`は同梱IPknotのnative improved beamを継続する。以下は切戻し前の導入実験と固定バイナリの評価記録であり、現在の既定動作は[README](../README.md)を参照。

## 適用内容

DAFS本体の `wip/linearization`、変更前 `8b62396646aec07d4b7d065baad21a7de1a984e3` を基準に評価した。
IPknot dev `fb0071d6629b06027d936130c8d724a740826e40` の再利用可能デコーダーと共有DP・ソルバー依存ソースを `vendor/ipknot` にコピーした。
上流の31ファイルには編集を加えておらず、SHA-256は `vendor/ipknot/UPSTREAM.json`、ライセンスは同ディレクトリの `COPYING` に記録している。
CMakeはこの同梱ソースをデフォルトとし、別のIPknot cloneや `IPKNOT_SOURCE_DIR` の指定を必要としない。

通常の `--fold-decoder=Nussinov` がこれらの実装を使用する。このオプションもデフォルトなので省略できる。
単一レベルの `DDNussinov` オラクルを呼ぶため、ILPソルバーもIPknot内側の多段DD反復も必要ない。
`--fold-decoder=IPknot` は別の多段・擬似結び目デコードを選ぶオプションであり、今回の切替には使っていない。

非線形の `SparseNussinov` は基準版8b62396でもすでにIPknotの厳密DPを使用していたので、そのアダプターを保った。
beam=0でclosed-candidate sparsificationを使い、DAFSの最小塩基対幅、目的関数、候補サポート、孤立塩基対方針を保持する。
線形の入力サイズによる旧DAFSカーネルへの切替を除去し、すべてIPknotのnative improved beamを使用する。
内部区間の全状態がbeamから落ちた場合には、上流どおり空の内部区間を許す。線形結果は変わり得る。
上流のpruning observerでDAFSの既存DD完了上界を計算し、停止判定の証明を維持する。

## 置換の診断

- 同一のランダム入力、負の乗数、同点スコア、孤立塩基対の許可・禁止、密・疎行列、最終プロファイルボーナスを含む1,848ケースで、非線形は変更前8b62396とスコアのfloatビットおよび構造が完全一致した。
- 線形264ケースのうち、構造またはスコアが変わったのは28ケース。目的関数は18ケース改善、7ケース低下、239ケース同点だった。有限beamには単調な改善保証はない。
- 独立した厳密DPと比較する最適値テスト、構造制約、負のスタック候補、DD上界、グラフ再利用、長い疎入力、serial/parallel一致のテストを通過した。
- OpenMP有効ビルドの全テストは30成功、ソルバーを無効化しているためILP専用1テストをskip。OpenMP無効のデコーダー3テスト、GLPK有効ビルドの追加2テストも成功した。
- 本体リポジトリを直接ソースにしたビルドも、外部IPknotパスなしで成功し、非線形・線形の通常CLI実行を確認した。

「完全一致」の比較基準は今回変更前の8b62396である。旧c098e3aや独立したdense参照実装に対しては、8b62396に含まれる既存のsparsification導入ですでに同点解や丸めに差が生じ得る。
今回、非線形に新たな同点解選択や丸め変更を加えていない。

## Murlet dataset1

85入力、旧推奨のN0-13/N1-34/L0-13/L1-15、MAFFT方式のDAFS距離ガイド木、孤立塩基対許可、max-iter=500、refinement=0。
各入力・条件を3回計測した。3回は時間計測の繰返しで、iterative refinementは無効。
TSUBAME課金グループtga-satolab-gtex、job 8997754、同一ノードr17n5、16並列、OMP_NUM_THREADS=1。
反復単位の前後実行順を before/after、after/before、before/after と交互にした。
SCIはViennaRNA 2.7.2、SPS/F1は既存ベンチマークの同一スコアラーを使用した。

|群|版|SCI|SPS|F1|平均wall秒/入力|
|---|---|---:|---:|---:|---:|
|N0|変更前8b62396|0.769698|0.782612|0.683647|2.7511|
|N0|同梱IPknotデフォルト|0.769698|0.782612|0.683647|2.7294|
|N1|変更前8b62396|0.818118|0.750926|0.672016|3.9337|
|N1|同梱IPknotデフォルト|0.818118|0.750926|0.672016|3.9064|
|L0|変更前8b62396|0.785907|0.775591|0.680838|6.5272|
|L0|同梱IPknotデフォルト|0.785907|0.775591|0.680838|6.6466|
|L1|変更前8b62396|0.812116|0.768919|0.689319|10.1595|
|L1|同梱IPknotデフォルト|0.812116|0.768919|0.689319|10.2069|

速度比は変更前÷変更後で、1より大きいと変更後が高速。95%区間は入力単位のpaired bootstrapで算出した。

|群|速度比|95%区間|予測が変わった入力|
|---|---:|---|---:|
|N0|1.0080×|1.0069–1.0091|0/85|
|N1|1.0070×|1.0062–1.0078|0/85|
|L0|0.9820×|0.9779–0.9852|0/85|
|L1|0.9954×|0.9922–0.9981|0/85|

全2,040実行が成功。前後のガイド距離行列・順序付き木は1,020組すべて一致。非線形の出力は510組すべてバイト一致した。各版内の反復間の予測もすべて一致した。

今回の4設定では精度向上を観測していない。IPknotのカーネルを使うことだけで、確率計算や外側DDを含む全体時間の改善は保証されない。表の速度はデコーダー単体ではなく予測器全体のwall時間である。

## 再現記録

- [集計JSON](ipknot_nussinov_bundled_default_2026-10-13.json)
- [集計CSV](ipknot_nussinov_bundled_default_2026-10-13.csv)
- 実験設定・全結果SQLite・実行ログ・診断用C++・集計スクリプトは、同じ親ディレクトリのdafs-benchmarkリポジトリにある `build/dafs-ipknot-nussinov-default/` に保存した。
- 変更前のソースはその下の `baseline-source/` に8b62396をarchiveして保存した。比較実行ファイルのSHA-256は集計JSONに記録している。
- 本体ソースからビルドした実行ファイルは `dafs-benchmark/build/dafs-ipknot-nussinov-default/native-cmake/src/dafs`。

デフォルトを使うビルド例（ViennaRNAのpkg-configパスは環境に応じて設定）:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_ILP=OFF
cmake --build build -j
```

既存CMakeキャッシュに外部IPKNOT_SOURCE_DIRが残っている場合は、新規buildディレクトリを使うか `cmake -U IPKNOT_SOURCE_DIR` でキャッシュ値を消す。
