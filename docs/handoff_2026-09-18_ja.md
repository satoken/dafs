# DAFS 引き継ぎメモ

作成日: 2026-09-18
対象ブランチ: `wip/linearization`
現在の HEAD: `709b5c9` (`origin/wip/linearization` と一致)

## 1. このメモの範囲

この作業ツリーは 2025-12-01 に取得され、`wip/linearization` の開始点は
`7757b6b`（dynamic CBP generation の実装）だった。以下は、その開始点から
現在の HEAD までにこの計算機で加えた内容である。

対象コミットは次の 9 件。

```text
729d81a 2025-12-05 Update CMake configuration and enhance IP solver support
946500d 2025-12-05 Refactor dynamic CBP logging
a1b81fb 2025-12-07 Refactor dynamic CBP generation
e0f65c5 2026-07-19 Improve dynamic CBP separation and stale-column pruning
1f61412 2026-07-20 Linearize DAFS probability and decoding pipeline
a1e8514 2026-07-20 Strengthen linear bounds and score pair matches
c966081 2026-08-02 Add reproducible benchmark instrumentation
9360ad8 2026-08-03 Tighten and optimize certified linear bounds
709b5c9 2026-09-18 Extend linear DAFS scoring and benchmark suite
```

`7757b6b..HEAD` の差分は 84 ファイル、約 20,075 行追加・729 行削除。

## 2. 実装上の主な進展

### IP ソルバーとビルド

- CMake の最低バージョンを 3.18 に更新。
- IP バックエンドを CMake オプションで切り替え可能にした。
  `GLPK` がデフォルトで、ほかに `CPLEX`、`Gurobi`、`SCIP`、`HiGHS` を追加。
- `cmake/FindHiGHS.cmake` を追加し、`src/ip.cpp` を SCIP/HiGHS を含む
  バックエンド対応に整理。IP 変数の上下限も修正。
- CTest に線形モデル、線形 Nussinov、RIBOSUM、ベンチマーク、入力検証の
  テストを追加。

### Dynamic CBP

- dynamic CBP 生成の初期化・不要な経路を整理。
- `std::cout` ベースのログを `spdlog` に統一。
- CBP の分離、正の reduced-cost 列の追加、投影配列の重複除去を整理。
- 一定期間使われず、かつ RIBOSUM 寄与もない stale CBP を削除する仕組みを追加。
- 生成した列を含む線形デコーダ用のサポート集合を管理するようにした。

### 線形化パイプライン

- `LinearAlign` に `LinearTurboFold`、`CONTRAlign`、`ProbConsRNA` のスコアモデルを
  持たせ、線形確率計算と線形 DD デコーダを接続。
- `LinearNeedlemanWunsch`、疎行列、`LinearNussinov` を追加。
- 線形 fold/align の beam、certified upper bound、feasible repair を統合。
- 線形構成では RNAalifold のプロファイル計算を無効化し、平均 BPP と
  LinearPartition 系の処理を使う。
- 主なオプションは `--align-model LinearAlign`、`-s lpc`、
  `--align-beam`、`--align-dd-beam`、`--linfold-beam`、
  `--fold-dd-beam`、`--fold-final-beam`、`--dynamic-cbp`。

### Certified bound と RIBOSUM

- `src/relaxed_bounds.h` に線形構造・アラインメントの certified bound を追加。
- beam bound、pruned-state bound、加法的 bound、row/column bound 等を DD trace に出力。
- 線形 Nussinov は pruning 状態を記録し、bound と feasible repair を使って停止判定する。
- `src/ribosum.{h,cpp}` を追加し、profile 間の pair-pair RIBOSUM85-60 スコアを実装。
- `--ribosum-weight` と、最終 Nussinov に使う `--final-ribosum-weight` を追加。
- `--max-iter=0` は IP バックエンドが必要であることを検証するようにした。

### ベンチマーク基盤

- `benchmarks/run.py` を追加・拡張。再開可能な per-run 結果、SHA-256、GNU time
  の wall/user/system time、最大 RSS、DAFS 内部 JSONL trace、環境 manifest を保存する。
- `benchmarks/score.py` で SPS、alignment PPV、structure sensitivity/PPV/MCC、
  pair-pair CBP precision/recall/F1、SCI を計算。
- ViennaRNA のバージョンを明示して SCI を記録するようにした。
- Murlet の manifest（872 件）と LinearTurboFold 非 viral データの manifest（65 件）、
  入力生成・検証スクリプトを追加。
- factorial、16S/23S scaling、RIBOSUM weight sweep、LinearTurboFold 比較用の
  JSON 設定、Slurm/TSUBAME の submit/run/merge スクリプトを追加。
- LinearTurboFold は DAFS の RIBOSUM 対応とは別の baseline として評価する。

### 文書

- [ベンチマーク手順](../benchmarks/README.md)
- [データ準備手順](../benchmarks/data/README.md)
- [線形化の日本語 draft](linearization_paper_draft_ja.md)
- [線形化説明スライド HTML](dafs_linearization_slides.html)
- `README.md` に線形モデル、RIBOSUM、最終 Nussinov bonus の説明を追記。

## 3. 保存済みベンチマーク結果

結果は Git 管理外（`.gitignore` の `benchmarks/results*/`）だが、この作業ツリーには残っている。
`analysis.json` の集計は次のとおり。

| 結果ディレクトリ | 結果 | 実行数 | Slurm job |
|---|---:|---:|---:|
| `results-comparison/dafs-accuracy` | 成功 | 156/156 | 296963 |
| `results-comparison/dafs-16s` | 成功 | 100/100 | 296964 |
| `results-comparison/dafs-23s` | 成功 | 20/20 | 296965 |
| `results-comparison/linearturbofold-16s` | 成功 | 25/25 | 296967 |
| `results-comparison/linearturbofold-23s` | 成功 | 5/5 | 296968 |
| `results-linear-ribosum-03/16s` | 成功 | 25/25 | 296970 |
| `results-linear-ribosum-03/23s` | 成功 | 5/5 | 296971 |
| `results-linear-contralign-ribosum-sweep-murlet` | 成功 | 195/195 | 297032/297033 |
| `results-linear-contralign-ribosum-sweep-16s` | 成功 | 375/375 | 297034/297063 |
| `results-linear-contralign-ribosum-sweep-23s` | 成功 | 75/75 | 297064/297065 |

未解決の結果は 2 件。

1. `results-comparison/linearturbofold-accuracy` は 39 件すべてプログラム実行自体は成功したが、
   scorer が `prediction and reference must share at least two sequences` で失敗した。
   LinearTurboFold 出力と Murlet reference の ID/配列対応を確認してから `--score-only` で再評価する。
2. `results-ribosum-sweep-murlet` は 390 件中 389 件が有効。
   `linear_ribosum_w0750 / murlet-dataset1-RF00003-4 / rep-01` の internal metrics が無効で、
   Slurm job 296969 は失敗終了になっている。該当 run の `stderr.log`、`metrics.jsonl`、
   `result.json` を確認し、必要なら `benchmarks/run.py --rerun-failed` で再実行する。

`results-rf00005` には 2026-08-01 の旧 smoke 実行結果がある。これは最新の factorial 集計とは別扱い。

## 4. ビルド・テスト状態

### 必須環境

- C++17、CMake 3.18 以上
- ViennaRNA の開発ファイル。`pkg-config --modversion RNAlib2` が成功すること
- IP を使う場合は GLPK（デフォルト）または CPLEX/Gurobi/SCIP/HiGHS
- ベンチマークは Python 3 と GNU `time` が必要
- Slurm 実験では `RNAalifold`、`RNAfold`、DAFS の x86_64 ビルドを compute node から参照できること

### 再ビルド例

```sh
VIENNA_PREFIX=/path/to/ViennaRNA
export PATH="$VIENNA_PREFIX/bin:$PATH"
export PKG_CONFIG_PATH="$VIENNA_PREFIX/lib/pkgconfig:$VIENNA_PREFIX/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="$VIENNA_PREFIX/lib:$VIENNA_PREFIX/lib64:${LD_LIBRARY_PATH:-}"
pkg-config --modversion RNAlib2

cmake -S . -B build-amd64/dafs \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$VIENNA_PREFIX" \
  -DBUILD_TESTING=ON
cmake --build build-amd64/dafs --parallel 8
ctest --test-dir build-amd64/dafs --output-on-failure
```

GLPK 以外を使う場合は `-DENABLE_CPLEX=ON`、`-DENABLE_GUROBI=ON`、
`-DENABLE_SCIP=ON`、`-DENABLE_HIGHS=ON` のいずれか一つを指定する。

### この計算機での確認結果

`/tmp/dafs-handoff-build` で CMake/CTest を実行したところ、GLPK と Python は見つかったが
ViennaRNA の `RNAlib2` が見つからず、DAFS 本体のビルドが `ViennaRNA/fold.h` で停止した。
そのため追加テスト 9 件のうち、依存しない次の 5 件だけが通過した。

```text
linear_models_fail_fast
relaxed_bounds
linear_nussinov
ribosum_pair_match
benchmark_scorer
```

DAFS 本体が生成されなかったため、benchmark runner smoke、negative max-iter、linear pipeline、
zero-max-iter の 4 件は実行不能。次の計算機では ViennaRNA の `pkg-config` 設定後に再実行する。

## 5. 次の計算機へ移すときの注意

1. `git fetch origin` 後、`wip/linearization` を checkout し、HEAD が `709b5c9` になることを確認する。
2. `build-amd64/dafs/src/dafs` は 2026-09-02 に作られた x86_64 バイナリで、現在の HEAD より古い可能性がある。
   成果物として信用せず、必ず再ビルドする。
3. Git 管理外の benchmark results、ダウンロード済みデータ、入力・reference は Git pull では移らない。
   必要なら `benchmarks/results-*` と `benchmarks/data` の ignored files を別途 rsync する。
4. 未追跡の資料成果物（スライド PNG/PDF/PPTX、`linearization_status_2026-07-23.pptx`）も必要なら別途コピーする。
   ソースとして再生成できるものは `docs/dafs_linearization_slides.html` と
   `docs/make_dafs_linearization_pptx.py` を使う。
5. データがない場合は `python3 benchmarks/data/prepare.py` を実行する。取得元、checksum、
   データの扱いは [データ準備手順](../benchmarks/data/README.md) を参照する。
6. Slurm 実験を再開する場合は、ViennaRNA のバージョン、PATH、GNU `time`、compute node の
   アーキテクチャを先に確認する。完了済み run は通常スキップされる。

## 6. 現在の作業ツリーの未追跡物

追跡対象の変更はない。未追跡なのは主に次のもの。

- `build-amd64/dafs/src/dafs` と CMake/CTest 生成物
- `docs/dafs_linearization_render/*.png`
- `docs/dafs_linearization_slides.{pdf,pptx}`
- `docs/dafs_linearization_slides_editable.pptx`
- `docs/dafs_linearization_slides_imagebacked_legacy.pptx`
- `linearization_status_2026-07-23.pptx`

これらを Git に入れるか、別途保存するかは次の計算機で判断する。特に build binary は再生成可能である。
