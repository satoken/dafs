# プロファイル熱力学項の RNAalifold 互換化と Murlet 比較

## 結論

LinearAlifold のプロファイル熱力学計算に、ViennaRNA の通常の特殊
tri-/tetra-/hexaloop 項を有効にするモードを追加した。Murlet dataset1 の
85ケースで比較したところ、SCI は4組すべてで改善した。一方、MCC は
`lpv × LinearAlign` で最も大きく改善したが、`lpc × LinearAlign` では低下した。
したがって、特殊 loop 項だけでMCC低下を一律に解消できるわけではない。

## 実装

`ProfileNearestNeighbor` に `EnergyMode::Legacy` と
`EnergyMode::RNAalifold` を追加した。DAFS のCLIでは次で切り替えられる。

```text
--linear-profile-energy=legacy
--linear-profile-energy=rnaalifold
```

RNAalifoldモードでは、次を行う。

- ViennaRNA 2.7.2 の37℃、`dangles=2` のモデルを使う。
- `special_hp=1` とし、各行のギャップを除いた配列から、短いヘアピンに
  closing pairを含む `size+2` の実配列を渡す。
- 既存の `vrna_E_hairpin`、`vrna_E_internal`、
  `vrna_E_multibranch_stem`、`vrna_E_exterior_stem` による評価を維持する。
- プロファイルの分配関数・outside・BPPはLinearFoldが計算し、DAFSの
  RIBOSUM85-60はこのPFには加えない。

旧モードは残して比較可能にした。CLIの既定値は `rnaalifold`、直接APIの
既定値は既存テストとの互換性のため `Legacy` である。

これは `vrna_pf_alifold()` をそのまま呼ぶ変更ではない。LinearFoldのbeam探索、
プロファイルの候補制約、ギャップを含む投影、outside計算はLinearFold側に残るため、
「RNAalifold互換」は熱力学プリミティブと特殊loop処理についての互換を意味する。
RNAalifoldの比較PF全体とビット単位で同一になることは主張しない。

## 検証

- `cmake --build build-arm64 -j2`: 成功
- `ctest --test-dir build-arm64 --output-on-failure`: 17/17成功
- 単一行プロフィール130ケースで、特殊loop表を有効にした独立ViennaRNA
  PF/BPPと `logZ` および全BPPを比較し、テスト許容差内で一致
- Murlet dataset1: 85ケース × 4組合せ × 2モード = 680 run
- 全680 runが成功し、`metrics_valid=false` と不完全評価は0件

ベンチマーク設定は
[`murlet-dataset1-profile-energy-comparison-arm64-v272.json`](../benchmarks/configs/murlet-dataset1-profile-energy-comparison-arm64-v272.json)、
全結果は
[`summary.csv`](../benchmarks/results-murlet-dataset1-profile-energy-comparison-arm64-v272/summary.csv)
に保存した。実行バイナリのSHA256は
`ad60c39ef6def8bd51aa0fb7d92649137cac2aba067bd3d159269160190ccaf6`。

条件は全て `max-iter=80`、全beam=100、`OMP_NUM_THREADS=1`、
`--alifold`（progressive/final）で固定した。数値は85ケースのmacro平均。

| 組合せ | profile energy | MCC | SCI | SPS | wall秒 | RSS KB |
|---|---|---:|---:|---:|---:|---:|
| lpc × LinearAlign | legacy | 0.683935 | 0.668361 | 0.746652 | 3.904 | 21822 |
| lpc × LinearAlign | RNAalifold | 0.681526 | 0.674943 | 0.746952 | 3.901 | 21798 |
| lpc × LinearCONTRAlign | legacy | 0.682858 | 0.679505 | 0.740201 | 3.884 | 22008 |
| lpc × LinearCONTRAlign | RNAalifold | 0.683348 | 0.684958 | 0.739714 | 3.885 | 22012 |
| lpv × LinearAlign | legacy | 0.669081 | 0.666066 | 0.760281 | 3.991 | 22697 |
| lpv × LinearAlign | RNAalifold | 0.677967 | 0.674773 | 0.761289 | 3.996 | 22671 |
| lpv × LinearProbCons | legacy | 0.671514 | 0.677710 | 0.748100 | 3.914 | 22068 |
| lpv × LinearProbCons | RNAalifold | 0.671479 | 0.679023 | 0.748566 | 3.918 | 22062 |

## legacy からの差

| 組合せ | ΔMCC | ΔSCI | ΔSPS | Δwall秒 | 構造が変化したケース |
|---|---:|---:|---:|---:|---:|
| lpc × LinearAlign | -0.002409 | +0.006583 | +0.000300 | -0.003 | 42/85 |
| lpc × LinearCONTRAlign | +0.000489 | +0.005453 | -0.000487 | +0.001 | 32/85 |
| lpv × LinearAlign | +0.008887 | +0.008707 | +0.001008 | +0.005 | 47/85 |
| lpv × LinearProbCons | -0.000035 | +0.001313 | +0.000466 | +0.003 | 38/85 |

RNAalifoldモードで予測構造が変わったケースは40〜52/85あり、差は丸め誤差だけではない。
ただし実行時間差は最大でも約5 ms/ケース、RSS差は約26 KBで、熱力学項の追加コストは
このデータセットでは実質的に無視できる。

## 解釈

今回の変更は、既存のプロファイル計算を別の熱力学モデルへ置き換えたというより、
既に使っていたViennaRNAのループ評価をRNAalifoldの通常設定に近づけ、特に
配列依存の特殊hairpin項を復活させたものだ。従って、MCCの差は特殊loopだけでなく、
そのスコア差がbeam候補・DD更新・最終デコーダーの境界を越えたことによって生じる。

今回の条件では `lpv × LinearAlign` がMCC/SPS/SCIの3指標で同時に改善した。一方、
`lpc × LinearAlign` はSCI/SPSがわずかに改善してもMCCが低下したため、実運用の既定値を
組合せごとに一律変更するのではなく、MCCを主目的にする場合は組合せ別に選ぶべきである。

ViennaRNAのRNAalifoldはアラインメント全体を比較モデルとして扱い、熱力学項に加えて
共変異項を持つ。今回のLinearFold経路では、共変異項は既存のRNAalifold互換スコアを使い、
DAFSのRIBOSUMはプロファイルPFから分離している。この区別はViennaRNAの
[RNAalifoldチュートリアル](https://viennarna.readthedocs.io/en/latest/tutorial/RNAalifold.html)
および[RNAalifoldマニュアル](https://viennarna.readthedocs.io/en/latest/man/RNAalifold.html)
の比較モデルと照合する際に重要である。
