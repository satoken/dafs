# Murlet dataset1: linear workspace再利用の効果

`BeamAlign::max_alignment()`のlayer/ranking/key領域と、
`LinearNussinov`のbeam/candidate/certificate/traceback領域を呼出し間で再利用した。
beam=100、max-iter=100または1000、85例、線形4条件の計680 runを、
変更前の同一ベンチマークと比較した。

## 同一性

- 成功: 680/680、internal metrics valid: 680/680
- prediction.faのバイト一致: 680/680
- MCC/SCI/SPSを含むquality辞書の一致: 680/680
- 時間以外のrun summary一致: 680/680

## 計算時間とメモリ

wall/DDは1例あたりの平均秒数。減少率が正なら高速化。RSSはGNU timeの最大常駐メモリ平均。

| condition | max-iter | MCC | SCI | SPS | wall old→new s | wall減 | DD old→new s | DD減 | RSS変化 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `cl_linear_ath02` | 100 | 0.6852 | 0.6792 | 0.7553 | 4.59→4.44 | 3.5% | 2.77→2.60 | 6.0% | +4.5% |
| `cl_linear_ath02` | 1000 | 0.6913 | 0.7445 | 0.7559 | 24.72→22.78 | 7.9% | 22.77→20.86 | 8.4% | +4.8% |
| `lpc_contralign_ath02` | 100 | 0.6857 | 0.6623 | 0.7536 | 4.50→4.36 | 3.2% | 2.60→2.45 | 5.5% | +4.1% |
| `lpc_contralign_ath02` | 1000 | 0.6933 | 0.7184 | 0.7529 | 22.34→21.43 | 4.1% | 20.56→19.46 | 5.3% | +4.6% |
| `lpv_linear_fpct05` | 100 | 0.6694 | 0.6809 | 0.7728 | 4.94→4.62 | 6.5% | 3.04→2.68 | 11.7% | +3.8% |
| `lpv_linear_fpct05` | 1000 | 0.6747 | 0.7642 | 0.7686 | 26.06→24.25 | 6.9% | 24.00→22.31 | 7.0% | +4.3% |
| `lpv_probcons_fpct05` | 100 | 0.6635 | 0.6609 | 0.7573 | 4.96→4.40 | 11.3% | 3.02→2.51 | 16.9% | +3.4% |
| `lpv_probcons_fpct05` | 1000 | 0.6703 | 0.7258 | 0.7566 | 24.95→23.12 | 7.4% | 23.09→21.20 | 8.2% | +3.6% |

## 解釈

- 精度・探索経路を変えず、DD内部時間は全条件で短縮した。
- 短縮幅はmax-iter=100で5.5〜16.9%、max-iter=1000で5.3〜8.4%。
- バッファを解放せず保持するため、最大RSSは3.4〜4.8%増えた。
- したがって今回の変更は、少量の常駐メモリと引き換えに反復ごとのallocator/hash再確保を削減する最適化である。

## 再現ファイル

- [設定](../benchmarks/configs/murlet-dataset1-buffer-reuse-arm64-v272.json)
- [比較CSV](../benchmarks/results-murlet-dataset1-buffer-reuse-arm64-v272/buffer_reuse_comparison.csv)
- [全run summary](../benchmarks/results-murlet-dataset1-buffer-reuse-arm64-v272/summary.csv)
