# LinearAlifoldありDD：alignment証明上界の改善（2026-09-20）

## 結論

alignment上界が緩くて止まれなかった短い診断例は解消できた。一方、Murlet全体では「上界が締まる」効果の割に反復削減は小さく、毎反復のblock計算費用もあるため、これだけで大幅な全体高速化にはならなかった。100反復で従来500反復と同等の品質に到達する改善ではない。

上界改善のみの全85ケース×4構成では、同じmax-iterの対照に対し、比較した1,020予測ファイルがすべてバイト一致した。MCC/SCI/SPSの差も0。既存の探索軌道を変えず、証明できたmergeだけ早く打ち切る設計と整合する結果である。ただし任意入力の出力一致を保証するものではない。

| 構成（max-iter=500） | 通常DD反復 | unpruned+block32 | 反復削減 | 通常gap | 改善gap |
|---|---:|---:|---:|---:|---:|
| A：lpc × LinearAlign | 4309.7 | 4050.6 | 6.0% | 2.635% | 1.146% |
| B：lpc × LinearAlign-CONTRAlign | 4329.8 | 4202.6 | 2.9% | 4.089% | 0.929% |
| C：lpv × LinearAlign | 4384.7 | 4206.5 | 4.1% | 3.010% | 1.295% |
| D：lpv × LinearAlign-ProbConsRNA | 4350.7 | 4277.3 | 1.7% | 3.521% | 0.662% |

DD反復は85ケース平均、各ケース9 mergeの合計。gapは各mergeの相対gapの平均。枝刈りなし上界だけの反復削減は構成別に約1.1〜2.1%。

MCC・SCI・SPS・時間・違反数を含む全条件の数値は[集計表](dafs_dd_certificates_tables_2026-09-20_ja.md)を参照。違反数は最終DD反復のsubproblem間不一致であり、回収後の返却解の違反数ではない。上界が締まって早く止まると、この値が僅かに増えることもある。

次に費用対効果を検討するなら、block上界を毎反復でなく一定間隔に限定する方法が候補。任意の反復で得た有効上界は履歴の最小として保持できるため、評価を省略しても上界の正当性は失われず、停止検出が遅れるだけである。あるいはblock幅内のprefix更新を小さな固定幅木で短縮する余地がある。どちらも今回は未実装で、今回の速度改善としては数えない。

## 実装の位置づけ

前報 `dafs_dd_recovery_2026-09-20_ja.md` の次の優先事項を実装して評価した。主対象は、良い解・乗数が得られてもalignmentの軽量上界が緩く、証明停止できない問題である。新機能はすべて既定off。モデル、事後確率、RIBOSUM重み、beam幅、DDの停止閾値は変更していない。

| オプション | 変更内容 |
|---|---|
| `--dd-unpruned-bound` | BeamAlignが実際に捨てた状態数を記録し、ゼロだった呼び出しだけ全探索のalignment最適値として上界に利用 |
| `--dd-block-bound=32` | 固定幅32のrow block内で単調alignmentを厳密に解き、block間制約を除いた緩和の上界を利用。beam/certified両方の乗数で評価 |
| `--dd-outward-lb` | 返却候補の構造・alignment・CBP整合を確認し、元の係数を下向きに再採点して実行可能下界を更新 |

block幅は0=off、1〜64。線形alignment/構造デコーダー・疎な乗数・有限かつ非負のalignment閾値を要求する。outward-LBは非負の有限な構造重みも要求する。

## 上界の根拠と計算量

### 枝刈りなしの検出

`BeamAlign::max_alignment()` は対角layerごとの単調alignment DPであり、別のenvelope制限はない。この呼び出し中にbeam上限で状態を一度も捨てていなければ、完全なDPである。単に「beam幅が大きそう」「tracebackに枝刈りがない」とは判定しない。各呼び出しでカウンターを初期化し、回収処理等でデコーダーを再利用する前に結果を保持する。

追加DPは行わない。既存の固定beam探索を維持し、追加コストは枝刈り時のカウンター更新と値の取得程度。scoreはfloatに縮める前のdoubleを保持し、非負のmatch和のdouble丸め誤差を外向きに補った後、float上界へ変換する。

### 固定幅block

任意の全体単調alignmentを各row blockへ制限すると、それぞれblock内で実行可能である。したがって各blockの最適値の和は全体最適値以上になる。blockをまたぐ列の重複・順序は許容するが、block内では同じ列・同じrowの再使用を禁止する。同じ列の更新は凍結したprefixから行う。

支持集合は確率支持と両trackの非ゼロ乗数支持の和集合。確率がゼロでも乗数だけが存在する座標を落とさない。重複はmax演算で吸収する。列番号を4回のradix passで並べるため、`O(M log L)` のソートや大域Fenwick木は使わない。幅をW、支持サイズをMとすると時間 `O(W M + L)`、メモリ `O(M + L + W)`。Wと疎な支持次数を固定した従来の前提下で配列長の線形性を維持する。

既存のrow/column上界とblock上界の小さい方を採用する。certified trackの勾配は従来のrow緩和解のままで、block解を異なる目的関数のsubgradientとして流用しない。動的CBP pricingと乗数射影も変更していない。全体の上界組立てには、丸め済み総和からのfloat減算を避け、doubleで上向き加算する経路を追加した。

## 短い固定問題の検査

元の6組×4構成に対し、上界改善の5方式を120実行、下界改善の4方式を96実行。長さ8〜9の同じ固定問題について、整数最適値の列挙と左端/row緩和LPを用い、前報と同じ1e-4の数値許容で上下界を検査した。これは浮動小数点演算全体の形式検証ではない。

| 固定問題 | 従来の反復数 | 改善後 | 効いた変更 |
|---|---:|---:|---|
| lpc × LinearAlign / mixed | 500 | 93 | 枝刈りなし上界 |
| lpv × LinearAlign / mixed | 500 | 36 | 枝刈りなし上界 |
| lpv × LinearAlign-ProbConsRNA / mixed | 500 | 63 | 枝刈りなし上界 |
| lpc × LinearAlign-CONTRAlign / stem_length | 500 | 185 | block32上界をcertified乗数にも適用 |
| lpc × LinearAlign / stem_shift | 500 | 1 | 下向き再採点 |

mixedは `GCGCGCGCG` / `CGCGCGCGC`。stem_lengthは `GGGAAACCC` / `GGAAACCC`。stem_shiftは `GCGAAACGC` / `GGGAAACCC`。500は収束せず反復上限へ到達したことを表す。改善後は既存のcertified-gap条件で停止した。

stem_shiftの整数最適値は約12.8310113726。従来はLB=12.8310127258がこれを上回り、UB=12.8310117722とのgapが−9.54e−7となっていた。下向き再採点ではLB=12.8310108185、gap=+9.54e−7となる。同じ予測のまま停止し、負gapを0へ切り上げる処理も、停止許容幅の緩和も行っていない。

ただし、下界修正だけで全上界の機械精度の問題が解決したとはいえない。今回の短い診断にも、従来の上界が整数最適値を約2.14e−7下回る例が残った。既存のfloat係数生成・CBP reduced-costの判定・構造上界も含む全経路の区間演算化は未実施である。今回追加した下界は、保存された有限float係数を実数と解釈した目的に対し保守的に再採点する。下界をPolyak更新にも渡すため、ごく小さい差でも以降の探索・同点選択が変わりうる。

さらに `dd_rational_certificate_check.py` で、保存float係数を有理数へ正確に変換し、24固定問題の整数最適値を誤差なしで列挙した。下界追加3方式×24問題の72下界はすべて許容誤差ゼロで最適値以下。対して従来の下界は24問題中18で最適値を僅かに上回った。上界についてはnormal/outward単独で各1問題が僅かに最適値を下回り、unpruned_outward/hybrid32_outwardのこの24問題では違反がなかった。有限なテスト集合の結果を任意入力の完全な浮動小数点証明とは解釈しない。

下界修正のMurlet予備384実行では、500反復時の48ケース中5ケースで品質指標が変わった（修正単独）。DのRF00017-1ではMCC −0.02110、SCI −0.11067、SPS +0.00161で、12ケース平均のSCI差も−0.00922となった。数値的な下界の正当化と、生物学的精度の向上は別である。この修正は既定offの独立オプションに留め、上界のみの効果と混同しない。

## Murlet評価条件と保存物

上界改善の予備評価は同じ12ケース×4構成×5方式×100/500反復=480実行。全体評価はdataset1の85ケース×4構成について、通常100/500、枝刈りなし500、枝刈りなし+block32の100/500を比較する。下界修正は別バイナリの同条件対照を伴う予備評価として分離した。

構成は前報のA=lpc×LinearAlign、B=lpc×LinearAlign-CONTRAlign、C=lpv×LinearAlign、D=lpv×LinearAlign-ProbConsRNA。全てalifoldあり、5種類のbeamは100、seed=42、refinement=0、eta=0.5、bp-updateなし。RIBOSUM等を含む完全な引数は保存configの各conditionに記録した。SCI評価にはViennaRNA 2.7.2を使用。

|構成|weight|DD RIBOSUM|最終RIBOSUM|align-pct|align-th|
|---|---:|---:|---:|---:|---:|
|A|4|0|0|0.25|0.005|
|B|2|0.3|0.005|0.5|0.01|
|C|4|0.01|0|0.5|0.01|
|D|2|0.3|0|0.25|0.01|

- 上界評価root：`benchmarks/results-dd-certificates-20260920/`
- 下界追加評価root：`benchmarks/results-dd-certificates-numeric-20260920/`
- 上界評価バイナリSHA256：`77efedf0047744c25fac88510b40fca91ef29f9803f8dd4318a85ce98c178607`
- 下界追加バイナリSHA256：`740a2201d0403e4a3043833936c23d20795399fb6fbeeaf4b50ad5a1d2ed80ea`
- 各root内に凍結binary、source hashes、config、予測、metrics、集計CSV/JSONを保存。
- 評価driver：`benchmarks/benchmark_dd_certificates.py`
- 新規テスト：`tests/dd_certificate_test.cpp`（2,000乱数問題×7幅、独立dense oracleと比較）、`tests/test_dd_certificates.py`（CLI、既定off、短い問題の停止・同一予測、負gap再現）。既存BeamAlignの差分テストにも枝刈り状態の初期化検査を追加。

上界評価と下界予備評価の一部は同時にjobs=16で実行したため、bulk wall時間には資源競合が含まれる。方式間の正確な速度比は、別途逐次実行した再計測と区別する。MCC/SCI/SPSと反復数はこの実行順序に依存しない。同じ85ケースを過去のパラメータ探索にも使用しており、独立検証集合での性能保証ではない。

逐次再計測の2ケース（RF00001-2 / RF00059-2）では、unpruned単独の全体時間は従来の約0.994〜1.005倍でほぼ同じ。hybrid32は約1.065〜1.105倍で遅かった。これらは枝刈りなし判定による反復削減がほぼなかったケースであり、上界計算自体の費用を見ている。短いRF00031-2を追加した別の逐次比較も集計表に掲載した。いずれも3回の中央値で、全85ケースの速度保証ではない。

短いRF00031-2では、unpruned単独の全体時間比がA=0.920、B=0.984、C=0.974、D=0.961で、約1.6〜8.0%短縮した。hybrid32はこのケースでunpruned以上の反復削減がなく、時間比1.100〜1.216と悪化した。枝刈りなしの判定は追加費用が小さく、適用できる長さ域で効果が出る。一方、短い配列で既に厳密alignment値が得られるときにblockを重ねて計算する費用も、この比較では含まれている。

CTestは17/17成功。上界機能を既定有効にする変更は行っていない。全体精度を維持して安価に試す候補はunpruned単独であり、毎反復のblock32計算と下界再採点は、目的と影響を区別した実験用オプションとして扱う。

block演算単体の固定支持次数8での検査では、幅32の中央値は長さ1,000で0.378ms、2,000で0.778ms、4,000で1.611ms、8,000で3.195ms、16,000で6.424ms。測定値も線形スケーリングと整合した。これはDAFS全体の時間測定とは別である。

再実行例：

```sh
python3 benchmarks/benchmark_dd_certificates.py prepare
python3 benchmarks/run.py --config benchmarks/results-dd-certificates-20260920/pilot-config.json --jobs 16
python3 benchmarks/benchmark_dd_certificates.py collect-pilot
python3 benchmarks/benchmark_dd_certificates.py diagnostics
python3 benchmarks/benchmark_dd_certificates.py full
python3 benchmarks/run.py --config benchmarks/results-dd-certificates-20260920/full-config.json --jobs 16
python3 benchmarks/benchmark_dd_certificates.py collect-full
python3 benchmarks/benchmark_dd_certificates.py audit
python3 benchmarks/benchmark_dd_certificates.py timing
python3 benchmarks/benchmark_dd_certificates.py timing-short
python3 benchmarks/benchmark_dd_certificates.py --numeric prepare
python3 benchmarks/run.py --config benchmarks/results-dd-certificates-numeric-20260920/pilot-config.json --jobs 16
python3 benchmarks/benchmark_dd_certificates.py --numeric collect-pilot
python3 benchmarks/benchmark_dd_certificates.py --numeric diagnostics
python3 benchmarks/dd_rational_certificate_check.py
python3 benchmarks/summarize_dd_certificates.py
```

prepare/diagnostics/full/timingは既存保存物を上書きしない。現在のsourceでprepareすると下界機能込みの新バイナリが凍結されるので、上記の旧SHAと一致する再解析には保存済みbinary/configを使用する。
