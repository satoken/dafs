# Murlet dataset1 新LinearAlifold max-iter sweep

新しいLinearAlifoldモジュールとRNAalifold互換プロファイル熱力学モードを固定し、`max-iter=100,200,...,1000`だけを変更した。
Murlet dataset1の85例、4条件×10段階、計3,400 runを評価した。全beamは100固定。

## 集計方法

- MCC、SCI、SPSは85例のmacro average。計算時間はrun単位のwall timeとDD内部時間の平均。
- actual DD iterはmetricsの`run_summary.dd_iterations`の平均。
- 制約違反数は各progressive mergeの最終DD反復に記録された`dd_summary.violated`。
  violation/caseは全mergeの違反数を1例単位で合計してから平均し、violation/mergeはmerge数で割った値。
- certified violation/mergeはlinear版のcertified dual trackにおける各merge最終反復値。

全runの状態: status={'success': 3400}, metrics_valid={'True': 3400}, quality_complete={'True': 3400}。

## 結果

### lpc × LinearCONTRAlign (`lpc_contralign`)

| max-iter | MCC | SCI | SPS | wall s/case | DD s/case | actual DD iter | violation/case | violation/merge | max violation | certified violation/merge |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.6835 | 0.6967 | 0.7386 | 5.75 | 5.73 | 899.6 | 203.56 | 22.62 | 729 | 32.45 |
| 200 | 0.6858 | 0.7348 | 0.7363 | 6.39 | 6.37 | 1791.6 | 156.80 | 17.42 | 517 | 26.72 |
| 300 | 0.6869 | 0.7401 | 0.7351 | 6.75 | 6.73 | 2671.1 | 138.31 | 15.37 | 520 | 23.66 |
| 400 | 0.6888 | 0.7472 | 0.7353 | 7.25 | 7.22 | 3540.2 | 113.64 | 12.63 | 427 | 22.85 |
| 500 | 0.6891 | 0.7524 | 0.7342 | 8.18 | 8.14 | 4395.4 | 105.89 | 11.77 | 169 | 21.15 |
| 600 | 0.6872 | 0.7545 | 0.7335 | 8.89 | 8.85 | 5239.3 | 102.71 | 11.41 | 226 | 20.61 |
| 700 | 0.6874 | 0.7583 | 0.7342 | 9.29 | 9.24 | 6078.1 | 95.27 | 10.59 | 248 | 20.44 |
| 800 | 0.6878 | 0.7597 | 0.7341 | 9.80 | 9.75 | 6913.2 | 93.09 | 10.34 | 126 | 19.61 |
| 900 | 0.6873 | 0.7628 | 0.7338 | 10.13 | 10.08 | 7747.3 | 88.66 | 9.85 | 123 | 19.33 |
| 1000 | 0.6876 | 0.7616 | 0.7334 | 11.05 | 10.99 | 8580.7 | 87.34 | 9.70 | 71 | 19.22 |

### lpc × LinearAlign (`lpc_linear`)

| max-iter | MCC | SCI | SPS | wall s/case | DD s/case | actual DD iter | violation/case | violation/merge | max violation | certified violation/merge |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.6884 | 0.6990 | 0.7449 | 5.79 | 5.77 | 897.8 | 257.07 | 28.56 | 621 | 40.00 |
| 200 | 0.6893 | 0.7317 | 0.7456 | 6.49 | 6.47 | 1783.9 | 191.76 | 21.31 | 502 | 32.28 |
| 300 | 0.6905 | 0.7482 | 0.7433 | 6.85 | 6.83 | 2656.4 | 167.52 | 18.61 | 455 | 30.11 |
| 400 | 0.6891 | 0.7534 | 0.7429 | 8.06 | 8.03 | 3519.4 | 145.19 | 16.13 | 381 | 29.25 |
| 500 | 0.6907 | 0.7578 | 0.7431 | 8.26 | 8.23 | 4368.3 | 137.15 | 15.24 | 366 | 28.07 |
| 600 | 0.6908 | 0.7607 | 0.7435 | 8.58 | 8.54 | 5207.8 | 128.96 | 14.33 | 365 | 27.10 |
| 700 | 0.6922 | 0.7632 | 0.7448 | 9.44 | 9.39 | 6043.0 | 120.26 | 13.36 | 346 | 26.60 |
| 800 | 0.6935 | 0.7701 | 0.7443 | 10.07 | 10.02 | 6871.4 | 109.06 | 12.12 | 311 | 25.35 |
| 900 | 0.6937 | 0.7716 | 0.7444 | 10.61 | 10.56 | 7698.4 | 104.81 | 11.65 | 276 | 25.98 |
| 1000 | 0.6942 | 0.7764 | 0.7440 | 11.36 | 11.31 | 8521.8 | 94.45 | 10.49 | 98 | 25.36 |

### lpv × LinearAlign (`lpv_linear`)

| max-iter | MCC | SCI | SPS | wall s/case | DD s/case | actual DD iter | violation/case | violation/merge | max violation | certified violation/merge |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.6737 | 0.6906 | 0.7561 | 5.94 | 5.92 | 899.1 | 338.18 | 37.58 | 593 | 52.27 |
| 200 | 0.6697 | 0.7316 | 0.7503 | 6.24 | 6.22 | 1786.3 | 276.93 | 30.77 | 756 | 41.61 |
| 300 | 0.6679 | 0.7500 | 0.7464 | 7.69 | 7.66 | 2660.4 | 209.14 | 23.24 | 572 | 37.79 |
| 400 | 0.6694 | 0.7636 | 0.7469 | 7.99 | 7.96 | 3526.4 | 198.80 | 22.09 | 475 | 38.82 |
| 500 | 0.6732 | 0.7706 | 0.7474 | 8.50 | 8.46 | 4384.5 | 170.93 | 18.99 | 448 | 35.79 |
| 600 | 0.6741 | 0.7842 | 0.7470 | 9.24 | 9.20 | 5234.3 | 162.51 | 18.06 | 491 | 34.85 |
| 700 | 0.6736 | 0.7920 | 0.7454 | 10.25 | 10.20 | 6085.4 | 146.09 | 16.23 | 491 | 33.19 |
| 800 | 0.6748 | 0.7970 | 0.7453 | 10.96 | 10.91 | 6928.3 | 129.14 | 14.35 | 214 | 31.74 |
| 900 | 0.6741 | 0.7941 | 0.7457 | 11.46 | 11.41 | 7761.2 | 129.67 | 14.41 | 457 | 31.67 |
| 1000 | 0.6738 | 0.8014 | 0.7449 | 12.57 | 12.51 | 8596.5 | 129.09 | 14.34 | 504 | 29.51 |

### lpv × LinearProbConsRNA (`lpv_probcons`)

| max-iter | MCC | SCI | SPS | wall s/case | DD s/case | actual DD iter | violation/case | violation/merge | max violation | certified violation/merge |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.6739 | 0.6951 | 0.7487 | 5.47 | 5.46 | 897.3 | 247.26 | 27.47 | 378 | 37.35 |
| 200 | 0.6769 | 0.7176 | 0.7433 | 6.38 | 6.36 | 1784.0 | 184.09 | 20.45 | 240 | 30.13 |
| 300 | 0.6773 | 0.7432 | 0.7450 | 7.03 | 7.00 | 2654.2 | 156.54 | 17.39 | 457 | 27.54 |
| 400 | 0.6783 | 0.7482 | 0.7448 | 7.89 | 7.86 | 3512.1 | 138.32 | 15.37 | 220 | 26.13 |
| 500 | 0.6796 | 0.7512 | 0.7445 | 8.50 | 8.46 | 4355.9 | 130.00 | 14.44 | 294 | 25.13 |
| 600 | 0.6818 | 0.7600 | 0.7448 | 9.33 | 9.29 | 5193.0 | 113.56 | 12.62 | 125 | 24.17 |
| 700 | 0.6817 | 0.7631 | 0.7446 | 10.08 | 10.04 | 6025.0 | 110.01 | 12.22 | 238 | 23.80 |
| 800 | 0.6814 | 0.7617 | 0.7444 | 10.48 | 10.43 | 6852.0 | 106.44 | 11.83 | 170 | 23.91 |
| 900 | 0.6815 | 0.7623 | 0.7443 | 11.26 | 11.21 | 7675.4 | 105.27 | 11.70 | 140 | 23.34 |
| 1000 | 0.6822 | 0.7663 | 0.7449 | 11.37 | 11.31 | 8497.4 | 96.14 | 10.68 | 124 | 22.90 |

## 固定した条件

- `lpc_contralign`: `-a LinearAlign-CONTRAlign -s lpc --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0.005 --align-pct=0.5 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100 --align-th=0.02 --linear-profile-energy=rnaalifold`（`--max-iter`のみ各行の値に置換）
- `lpc_linear`: `-a LinearAlign -s lpc --alifold --dynamic-cbp --ribosum-weight=0.5 --final-ribosum-weight=0.01 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100 --align-th=0.02 --linear-profile-energy=rnaalifold`（`--max-iter`のみ各行の値に置換）
- `lpv_linear`: `-a LinearAlign -s lpv --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --align-pct=0.5 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100 --fold-pct=0.5 --linear-profile-energy=rnaalifold`（`--max-iter`のみ各行の値に置換）
- `lpv_probcons`: `-a LinearAlign-ProbConsRNA -s lpv --alifold --dynamic-cbp --ribosum-weight=0.3 --final-ribosum-weight=0 --weight=4 --max-iter=100 --align-beam=100 --align-dd-beam=100 --linfold-beam=100 --fold-dd-beam=100 --fold-final-beam=100 --fold-pct=0.5 --linear-profile-energy=rnaalifold`（`--max-iter`のみ各行の値に置換）

## ファイル

- [再現用設定](../benchmarks/configs/murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272.json)
- [全runのsummary.csv](../benchmarks/results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272/summary.csv)
- [集計CSV](../benchmarks/results-murlet-dataset1-linearalifold-profile-max-iter-sweep-arm64-v272/max_iter_summary.csv)
