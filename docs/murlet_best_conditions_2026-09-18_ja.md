# Murlet 条件別の最高 MCC

Murlet dataset1 の 13 データセット、各 3 反復（39 run）の確認結果を
平均した。各組み合わせについて、条件探索で MCC が最も高かった条件を
選んでいる。SCI は ViennaRNA 2.7.2 で計算した。選択した run はすべて
`metrics_valid=True` である。LinearTurboFold も同じ13データセット・3反復で
評価した。

| 組み合わせ | 条件ID | SPS | alignment PPV | 構造 sensitivity | 構造 PPV | MCC | CBP F1 | SCI |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| ViennaRNA × ProbCons | `vp_alifold_mcc` | 0.8172 | 0.8205 | 0.6611 | 0.7912 | **0.7140** | 0.6110 | 0.7551 |
| CONTRAfold × CONTRAlign | `cc_noalifold_mcc` | 0.8067 | 0.8093 | 0.6823 | 0.7640 | **0.7162** | 0.6039 | 0.7231 |
| lpc × LinearAlign | `cl_linear_mcc` | 0.8091 | 0.8190 | 0.6863 | 0.7625 | **0.7188** | 0.6124 | 0.7442 |
| lpv × LinearProbCons | `lpv_probcons_base` | 0.8171 | 0.8197 | 0.7011 | 0.7494 | **0.7202** | 0.6156 | 0.7487 |
| lpc × LinearCONTRAlign | `lpc_contralign_alignpct05` | 0.7881 | 0.7904 | 0.6705 | 0.7719 | **0.7146** | 0.6027 | 0.6932 |
| lpv × LinearAlign | `lpv_linear_mcc_apct05` | 0.8210 | 0.8259 | 0.7033 | 0.7510 | **0.7231** | 0.6169 | 0.7344 |
| LinearTurboFold | `linearturbofold` | 0.8452 | 0.8346 | 0.6733 | 0.6467 | 0.6559 | 0.5509 | 0.6255 |

## 再現用パラメータ

共通設定は `--dynamic-cbp`、`--max-iter=80`。線形条件では
`--align-beam=100 --align-dd-beam=100 --linfold-beam=100
--fold-dd-beam=100 --fold-final-beam=100` を使用した。`--weight`、
`--align-pct`、`--align-th`、`--fold-pct`、`--fold-th` の未指定値は
DAFSのデフォルト値（それぞれ 4、0.25、0.01、0.25、0.2）である。

```text
# ViennaRNA × ProbCons
-a ProbCons -s Vienna --alifold --dynamic-cbp
--ribosum-weight=0.02 --final-ribosum-weight=0
--fold-th=0.3 --max-iter=80

# CONTRAfold × CONTRAlign
-a CONTRAlign -s CONTRAfold --no-alifold --dynamic-cbp
--ribosum-weight=0.1 --final-ribosum-weight=0
--fold-pct=0 --max-iter=80

# CONTRAfold × LinearAlign
-a LinearAlign -s lpc --no-alifold --dynamic-cbp
--ribosum-weight=0.5 --final-ribosum-weight=0.01
--align-th=0.005 --max-iter=80
--align-beam=100 --align-dd-beam=100 --linfold-beam=100
--fold-dd-beam=100 --fold-final-beam=100

# lpv × LinearAlign
-a LinearAlign -s lpv --no-alifold --dynamic-cbp
--ribosum-weight=0.3 --final-ribosum-weight=0
--align-pct=0.5 --max-iter=80
--align-beam=100 --align-dd-beam=100 --linfold-beam=100
--fold-dd-beam=100 --fold-final-beam=100

# lpv × LinearProbCons
-a LinearAlign-ProbConsRNA -s lpv --no-alifold --dynamic-cbp
--ribosum-weight=0.3 --final-ribosum-weight=0
--weight=4 --max-iter=80
--align-beam=100 --align-dd-beam=100 --linfold-beam=100
--fold-dd-beam=100 --fold-final-beam=100

# lpc × LinearCONTRAlign
-a LinearAlign-CONTRAlign -s lpc --no-alifold --dynamic-cbp
--ribosum-weight=0.3 --final-ribosum-weight=0.005
--align-pct=0.5 --max-iter=80
--align-beam=100 --align-dd-beam=100 --linfold-beam=100
--fold-dd-beam=100 --fold-final-beam=100

# LinearTurboFold
hmm_beam=100 cky_beam=100 iterations=3
min_helix_length=3 pk_iterations=1 threshold=0.3
```

## 使用した結果

- `results-three-combo-confirm-arm64-v272/summary.csv`
- `results-two-linear-confirm-arm64-v272/summary.csv`
- `results-comparison/linearturbofold-accuracy/benchmark/summary.csv`

## 今回のLinearAlign追加条件探索

lpv × LinearAlign は、基準条件に対して `--weight=2/8`、
`--align-pct=0/0.5`、`--fold-pct=0/0.5`、`--align-th=0.005/0.02`、
`--fold-th=0.1/0.3` を比較した（11条件、13データセット、各1反復）。
全143 runで `metrics_valid=True` だった。MCC最良は
`--align-pct=0.5`（0.7231）で、SCIは0.7344だった。

今回の lpv × LinearProbCons は、`weight=2` ではなく、確認実験で
MCC 0.7202 となったデフォルト `weight=4` (`lpv_probcons_base`) を採用した。
