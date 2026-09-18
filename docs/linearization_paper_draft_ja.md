# 厳密な双対上界を保持する線形時間 DAFS と塩基対置換スコアの導入

> **詳細ドラフト（改訂版、2026-08-04）**
> 本稿は `wip/linearization` に実装済みの内容を、DAFS 原論文の記法に合わせて可能な限り網羅的に記述した初稿である。DAFS 原論文に現れない計算機科学・数理最適化上の用語は、初出時に日本語で意味と役割を説明する。論文投稿時には、未完了の統計検定、データセットの正式な出典情報、比較手法、図、実測スケーリング、バージョン番号を追記し、実装史に属する細部を本文または補足資料から削る必要がある。

## 要旨

RNA の配列アラインメントと共通二次構造予測を同時に解く DAFS は、塩基対確率とアラインメント確率に基づく最大期待精度（maximum expected accuracy; MEA）問題を整数計画として定式化し、構造とアラインメントの整合性制約をラグランジュ緩和することで、二つの RNA folding 問題、一つの pairwise alignment 問題、および consensus base-pair（CBP）選択問題へ分解する。しかし、従来実装には二次以上の確率計算、密な確率行列とラグランジュ乗数、二次または三次の動的計画、ならびに多数の CBP 候補が残り、長鎖 RNA への適用を制限していた。また、単に beam search へ置換した目的値は元の部分問題の最大値ではないため、ラグランジュ双対値を元問題の厳密な上界として解釈できなくなる。

本研究では、DAFS の分解構造を保ったまま、確率計算、復号、乗数更新、CBP 管理、および上界計算を線形時間の処理へ置換する。塩基対確率には LinearPartition、マッチ確率には LinearAlign/BeamAlign を用い、閾値以上の非零要素だけを保存する疎行列を全段で保持する。folding と alignment のラグランジュ部分問題には、各計算位置で有望な状態を一定個だけ残す beam search に基づく最大得点復号器を用いる。動的 CBP には、現在解から新しい候補を作る双方向候補分離に加え、まだ追加されていない候補が目的値を改善するかを表す reduced cost を漏れなく調べる pricing と、現在有用な候補だけを保持する working-set 管理を導入する。beam search の解は実行可能解として下界回復に用いる一方、厳密上界には、枝刈りされた状態が将来得られる最大得点を頂点被覆に基づいて評価する folding certificate、および alignment の制約を弱めた行・列容量緩和を用いる。さらに、これらとは独立に、構造の左端点容量緩和と alignment の行容量緩和からなる第二の双対最適化系列を保持し、二種類の上界の最小値を採用する。乗数の更新幅には、現在の上界候補と実行可能下界の差、および subgradient の大きさから決める Polyak step を用いる。

従来目的関数では係数がゼロであった pair-pair match 変数に RIBOSUM85-60 の塩基対置換スコアを付与し、その項を分解後の CBP reduced cost、pricing、下界 repair、および整数計画にも一貫して含める。Murlet 13 データセットにおいて、RIBOSUM 重み 0.075 は線形構成の平均 SPS を 0.8220 から 0.8286、MCC を 0.6502 から 0.6627、CBP F1 を 0.5690 から 0.5828 へ改善した。16S および 23S 長鎖データでは、RNAalifold を用いない非線形構成に対して線形構成はそれぞれ 1.82 倍、4.61 倍高速で、最大常駐メモリの平均をそれぞれ約 20%、48%削減した。以上により、固定閾値、固定 beam 幅、固定双対反復数、および固定配列数の条件下で、厳密な primal–dual bound を維持する長さ方向の線形時間・線形メモリ DAFS を実現した。

## 1. はじめに

相同 RNA のアラインメントでは、一次配列の類似性だけでなく、保存された二次構造、特に compensatory substitution を考慮することが重要である。配列を先に固定してから共通構造を予測する逐次法は、初期アラインメントの誤りを構造推定が修正できない。一方、構造とアラインメントを同時最適化する厳密な Sankoff 型動的計画は一般に高い時間・空間計算量を要する。

DAFS（Dual decomposition for Aligning and Folding Simultaneously）は、構造とアラインメントの事後確率を用いる MEA 推定を整数計画化し、構造整合性を双対分解することで、既存の folding decoder と alignment decoder を再利用可能にした手法である [Sato et al., 2012]。原論文は、二つの Nussinov 型部分問題、一つの Needleman–Wunsch 型部分問題、および独立な CBP 変数へ分解し、subgradient 法によりラグランジュ乗数を更新した。しかし、分解後の各部分問題や確率計算が長さに対して線形でなければ、分解そのものは長鎖 RNA に対する線形スケーリングを保証しない。

近年、固定幅 beam search により partition function と塩基対確率を近似する LinearPartition、および長鎖 RNA の確率的アラインメントに beam search を用いる LinearAlign/BeamAlign が提案された。これらは probability semiring 上の log-sum-exp 計算を固定 beam 幅で枝刈りする。本研究では確率推定のみならず、DAFS の max-sum 復号、双対変数、CBP 列生成、下界回復、上界証明までを一体として線形化する。

本研究の主な貢献は次の通りである。

1. LinearPartition と LinearAlign による確率推定、および閾値疎行列を DAFS の progressive alignment 全体へ導入した。
2. DAFS の folding/alignment 部分問題と最終共通構造予測を、確率とラグランジュ項を直接スコアとする固定 beam の max-sum decoder へ置換した。
3. 動的 CBP の heuristic separation、厳密 pricing、および stale-column pruning を組み合わせ、制限された列集合を用いても元問題に対する上界性を失わない working-set 法を構成した。
4. beam 解と厳密双対上界を分離し、folding の beam-pruning certificate、alignment の matching relaxation、および独立な certified dual track によって厳密上界を線形時間で計算した。
5. 合意塩基対 repair と best primal/dual 値の保持により、単調な実行可能下界と相対 gap を得た。
6. 従来係数ゼロであった pair-pair match 変数に RIBOSUM85-60 スコアを導入し、精度を改善した。
7. 機械可読 JSONL、資源使用量、精度、チェックサム、実行環境、primal–dual invariant を保存する再現可能なベンチマーク基盤を整備した。

## 2. 従来手法 DAFS

### 2.1 記法

原論文に従い、二つの RNA 配列または progressive alignment 中の二つの profile をそれぞれ (a=a_1\ldots a_{L_1})、(b=b_1\ldots b_{L_2}) とする。ここで profile とは、既に整列された一群の配列を一つの対象として扱ったものであり、profile の位置は alignment column に対応する。profile の配列数を (N_1,N_2) とする。構造アラインメントは、(a) 上の二次構造 (x)、(b) 上の二次構造 (y)、および両者のアラインメント (z) からなる。

- (x_{ij}\in\{0,1\}): (a) の位置 (i,j) が塩基対を形成する。
- (y_{kl}\in\{0,1\}): (b) の位置 (k,l) が塩基対を形成する。
- (z_{ik}\in\{0,1\}): (a_i) と (b_k) が同じ alignment column に配置される。
- (w_{ijkl}\in\{0,1\}): 塩基対 ((i,j)) と ((k,l)) が pair-pair match、すなわち同一の consensus base pair を表す。

添字は (i<j)、(k<l) とする。許容構造集合を ({\cal S}_1,{\cal S}_2)、単調な一対一アラインメント集合を ({\cal A})、有効な CBP 四つ組集合を ({\cal C}) と書く。擬似結び目を許さない場合、({\cal S}) は各位置の次数が高々 1 で、二つの塩基対 ((i,j),(i',j')) が (i<i'<j<j') のように交差しない集合である。

塩基対事後確率を (p^x_{ij}=P(x_{ij}=1\mid a))、(p^y_{kl}=P(y_{kl}=1\mid b))、マッチ事後確率を (p^z_{ik}=P(z_{ik}=1\mid a,b)) とする。構造閾値を (\theta_s)、アラインメント閾値を (\theta_a)、構造項の相対重みを (\omega) とする。現実装の profile-size 補正を含む局所係数は

\[
c^x_{ij}=\omega\frac{2N_1}{N_1+N_2}(p^x_{ij}-\theta_s),\quad
c^y_{kl}=\omega\frac{2N_2}{N_1+N_2}(p^y_{kl}-\theta_s),
\]

\[
c^z_{ik}=p^z_{ik}-\theta_a
\]

である。単一配列同士では (N_1=N_2=1) なので、原 DAFS の (\omega(p-\theta_s)) に一致する。複数の IPknot level を用いる場合は level ごとの (\theta_s) が存在するが、本稿の線形 decoder と実験は単一閾値を主対象とする。

### 2.2 MEA 目的関数と整数計画

DAFS は二つの構造とアラインメントの事後分布を独立と近似して期待 gain を因数分解する。従来目的関数は

\[
\max_{x,y,z,w}
\left\{
\sum_{i<j}c^x_{ij}x_{ij}
+\sum_{k<l}c^y_{kl}y_{kl}
+\sum_{i,k}c^z_{ik}z_{ik}
\right\}
\tag{1}
\]

であり、原実装では (w_{ijkl}) 自体の係数はゼロであった。制約は (x\in{\cal S}_1)、(y\in{\cal S}_2)、(z\in{\cal A}) に加え、各 consensus base pair の構造・アラインメント整合性を課す。実装と同値な形で、

\[
x_{ij}=\sum_{(k,l):(i,j,k,l)\in{\cal C}}w_{ijkl},
\tag{2}
\]

\[
y_{kl}=\sum_{(i,j):(i,j,k,l)\in{\cal C}}w_{ijkl},
\tag{3}
\]

および各 match ((i,k)) について

\[
\sum_{(i,j,k,l)\in{\cal C}}w_{ijkl}
+\sum_{(j,i,l,k)\in{\cal C}}w_{jilk}
\le z_{ik}
\tag{4}
\]

を課す。すなわち、選択された塩基対は相手 profile のちょうど一つの塩基対と対応し、その両端がアラインされなければならない。構造制約により一塩基は高々一塩基対へ参加するため、式 (4) は endpoint ごとの容量制約となる。

### 2.3 ラグランジュ緩和と双対分解

式 (2)–(4) を緩和する。式 (2),(3) の等式乗数をそれぞれ (q^x_{ij},q^y_{kl}\in\mathbb{R})、式 (4) の不等式乗数を (q^z_{ik}\ge0) とする。本実装の符号規約ではラグランジアンは

\[
\begin{aligned}
L(x,y,z,w;q)=
&\sum_{i<j}(c^x_{ij}-q^x_{ij})x_{ij}
+\sum_{k<l}(c^y_{kl}-q^y_{kl})y_{kl}\\
&+\sum_{i,k}(c^z_{ik}+q^z_{ik})z_{ik}\\
&+\sum_{(i,j,k,l)\in{\cal C}}
(q^x_{ij}+q^y_{kl}-q^z_{ik}-q^z_{jl})w_{ijkl}.
\end{aligned}
\tag{5}
\]

固定した (q) に対する双対関数 (D(q)=\max L(x,y,z,w;q)) は、次の独立な部分問題へ分解する。

\[
D_x(q)=\max_{x\in{\cal S}_1}\sum_{i<j}(c^x_{ij}-q^x_{ij})x_{ij},
\tag{6}
\]

\[
D_y(q)=\max_{y\in{\cal S}_2}\sum_{k<l}(c^y_{kl}-q^y_{kl})y_{kl},
\tag{7}
\]

\[
D_z(q)=\max_{z\in{\cal A}}\sum_{i,k}(c^z_{ik}+q^z_{ik})z_{ik},
\tag{8}
\]

\[
D_w(q)=\sum_{(i,j,k,l)\in{\cal C}}
\max\{0,q^x_{ij}+q^y_{kl}-q^z_{ik}-q^z_{jl}\}.
\tag{9}
\]

したがって (D(q)=D_x+D_y+D_z+D_w) である。任意の双対実行可能 (q) に対し (D(q)) は元の最大化問題の上界であり、双対問題は (min_qD(q)) である。従来 DAFS は (6),(7) を Nussinov 型 DP、(8) を Needleman–Wunsch 型 DP、(9) を係数が正の (w) の独立選択で解き、構造と CBP の不一致を subgradient として乗数を更新した。

### 2.4 Progressive multiple alignment

原 DAFS は pairwise joint alignment/folding を UPGMA guide tree に沿って progressive に適用する。profile の塩基対確率は各配列の確率を現在の alignment column へ射影して平均し、profile 間の match probability も構成配列対の確率を平均して求める。各 internal node で二 profile を統合し、最後に得られた multiple alignment 上で共通二次構造を復号する。したがって、本研究で「長さに線形」と呼ぶ単位は各 profile merge および各配列の確率計算である。配列数 (M) 自体を増やす場合、全配列対の posterior や guide tree 構築には別途 (M^2) 依存が残り得る。

### 2.5 従来実装のボトルネック

従来構成では、partition function、alignment posterior、密な (L_1^2,L_2^2,L_1L_2) 行列、Nussinov/Needleman–Wunsch DP、および CBP の直積候補が線形ではない。特に CBP は概念上四添字 (w_{ijkl}) であり、無制限に列挙すれば最悪 (O(L_1^2L_2^2)) となる。また、beam decoder の返す値をそのまま (6)–(8) の最大値とみなすと、値は真の最大値以下であり、式 (5) の双対上界性が失われる。この問題は速度と bound certificate を同時に扱う必要がある。

### 2.6 本研究で新たに用いる用語

以下の用語は DAFS 原論文の中心的な記法には含まれないため、本稿では次の意味で使用する。

**疎行列（sparse matrix）**は、行列の全要素を保存せず、閾値以上または非零の要素とその座標だけを保存するデータ構造である。長さ \(L\) の配列に対する塩基対確率を通常の行列で保存すると \(L^2\) 要素が必要になるが、一位置あたりの保存候補数が定数なら疎行列の要素数は \(O(L)\) となる。

**support（非零座標集合）**は、疎行列で実際に保存されている座標の集合である。例えば \(p^x\) の support は、塩基対確率が閾値を通過した \((i,j)\) の集合である。

**beam search** は、動的計画の各段階で全ての途中状態を残す代わりに、評価値が高い状態を高々 \(B\) 個だけ残す近似探索である。\(B\) を **beam 幅**と呼ぶ。固定 beam 幅なら保持状態数が配列長とともに増えないため高速化できるが、除外された状態の中に最適解へ至る状態が含まれる可能性がある。

**枝刈り（pruning）**は、以後の探索対象から途中状態または候補を除外する操作である。beam 幅を超えたために行う枝刈りは近似を生む。一方、本稿の suffix dominance のように、別の状態が常に同等以上であることを証明して除く操作は最適解を失わない。

**log-sum-exp** は、log probability \(a_r\) から \(\log\sum_r\exp(a_r)\) を数値的に安定して計算する演算であり、partition function や事後確率の計算に用いる。**max-sum** は総和を最大値へ置き換え、最も得点の高い一つの解を求める演算である。「semiring」はこれらの演算規則の組を指す代数学上の呼称だが、本稿では原則として「確率の総和計算」と「最大得点計算」と記述する。

**復号器（decoder）**は、確率または局所得点から、制約を満たす一つの構造またはアラインメントを選ぶアルゴリズムである。本稿の folding decoder は非交差塩基対集合を、alignment decoder は単調な一対一マッチ集合を返す。

**working set（作業候補集合）**は、全候補のうち現在の反復で明示的に保持して最適化に使用する候補集合である。CBP 全体を \({\cal C}\)、反復 \(t\) の working set を \({\cal C}_t\) とする。

**候補分離（separation）**は、現在解が示す構造とアラインメントの不一致から、新たに必要そうな CBP を見つける操作である。これは高速な発見法だが、それだけでは全ての有益な欠落候補を発見した保証はない。

**reduced cost（被約費用）**は、現在のラグランジュ乗数の下で、ある変数を 0 から 1 に変えたときにラグランジアンがどれだけ増えるかを表す係数である。本稿の CBP \((i,j,k,l)\) では式 (18) の \(\bar c_{ijkl}\) がこれに当たる。最大化問題なので、正の reduced cost を持つ CBP は双対部分問題で選ぶ価値がある。

**pricing（価格付け）**は、working set の外に正の reduced cost を持つ変数がないかを探索し、見つけた変数を working set へ追加する操作である。本稿で **exact pricing** と呼ぶのは、正の reduced cost を持つ欠落 CBP を漏れなく発見することを証明できる探索である。「exact」は確率計算そのものが厳密であるという意味ではなく、現在定義された CBP 候補集合に対する探索漏れがないという意味である。

**stale column（長期間不活性な列）**は、一定反復の間、解にも上界計算にも寄与しなかった CBP 候補である。column は整数計画の一変数、ここでは一つの \(w_{ijkl}\) を意味する。stale-column pruning はそのような変数を working set から除く操作である。

**主問題実行可能解（primal feasible solution）**は、構造、アラインメントおよび CBP の全整合性制約を満たす解である。最大化問題では、その得点は最適値の下界となる。**双対上界（dual upper bound）**は、制約をラグランジュ緩和した問題から得る最適値の上界である。両者の差を **primal–dual gap** と呼ぶ。gap が小さいほど、現在の解と証明可能な最適値の範囲が狭い。

**repair（実行可能解への修復）**は、互いに一致しない複数の部分問題解から、全制約を満たす主問題実行可能解を作り直す操作である。本稿では両構造が合意した塩基対だけを残す intersection repair と、固定したアラインメント上で共通塩基対を再最適化する consensus repair を用いる。

**relaxation（緩和）**は、一部の制約を取り除いて元問題より広い解集合を作ることである。最大化問題では緩和問題の最大値が元問題の上界になる。row-capacity relaxation は「各行から高々一つ」という制約だけを残し、アラインメントの単調性などを外す。

**certificate（上界証明値）**は、近似探索が返した得点とは別に、真の部分問題最適値がある値を超えないことを計算可能な不等式で保証する情報である。本稿の folding certificate は、枝刈りされた状態から将来得られる最大得点を上から評価する。

**vertex-cover potential（頂点被覆ポテンシャル）**は、各配列位置 \(i\) に非負値 \(u_i\) を割り当て、全候補塩基対 \((i,j)\) の得点 \(s_{ij}\) に対して \(u_i+u_j\ge s_{ij}\) を満たす値である。一つの塩基は高々一塩基対にしか使えないため、任意の構造得点は \(\sum_i u_i\) 以下になる。これはグラフ理論の重み付き頂点被覆の双対的な考え方に相当する。

**dual track（双対最適化系列）**は、固有のラグランジュ乗数と更新履歴を持つ一つの双対最適化過程である。本実装は、beam 解を用いて良い主問題解を探す通常トラックと、厳密に解ける緩和上界を直接小さくする certified track の二つを並行して保持する。

**Polyak step** は、現在の双対目的値と既知の主問題下界との差を subgradient の二乗ノルムで割って更新幅を決める方法である。上界候補と下界が近づけば更新が小さくなり、subgradient が大きければ一座標あたりの更新を抑える。

**外向き丸め（outward rounding）**は、上界を低い精度の数値型へ変換するとき、丸め後の値が計算値を下回らない方向、すなわち正の無限大方向へ丸めることである。浮動小数点誤差によって上界保証を失わないために用いる。

**cache（計算結果の再利用）**は、以前計算した値または候補一覧を保存し、同じ計算を繰り返さない仕組みである。cache は数学的な計算結果を変えず、実行時間だけを減らす。

## 3. 提案手法

### 3.1 線形性の定義と基本方針

本稿では、確率閾値 (\theta_s,\theta_a>0)、各 beam 幅 (B_p^s,B_p^a,B_d^s,B_d^a,B_f^s)、双対反復上限 (T)、profile に含まれる配列数を定数としたとき、profile 長 (L=L_1+L_2) に対して時間・追加メモリが (O(L)) となることを「線形」と呼ぶ。beam 幅や (1/\theta) は大きな定数になり得るため、これはパラメータを固定した場合の線形性である。また、LinearPartition/LinearAlign 自体は近似確率計算である。従って、「枝刈り前の熱力学・統計モデルを厳密に計算すること」と、「近似計算後に得られた疎な確率行列で定義される DAFS 目的関数の最適値を上から保証すること」を区別する。本研究が保証するのは後者である。

処理は、(i) 線形近似 posterior、すなわち各塩基対または match が正解に含まれる事後確率の計算、(ii) 閾値以上の確率だけを残す疎化、(iii) beam search による部分問題解と乗数更新方向の計算、(iv) 制約を弱めた問題による厳密上界、(v) 部分問題解を整合させる repair による下界、(vi) 欠落した有益な CBP を全て追加する exact pricing、からなる。

### 3.2 線形時間の確率計算と疎行列

folding model に `lpc` を指定した場合、CONTRAfold パラメータを用いる LinearPartition-C を、`lpv` または `LinFold` では Vienna/Turner パラメータを用いる LinearPartition-V を呼び出す。LinearPartition は、RNA の左から右へ部分構造を展開し、各位置で高得点の途中状態だけを固定数残すことにより、partition function と塩基対事後確率を近似する。標準 beam 幅は (B_p^s=100) である。

alignment model `LinearAlign` は LinearTurboFold 由来の BeamAlign 実装を用いる。forward 計算は配列先頭から各途中状態へ到達する全経路の確率を、backward 計算はその状態から配列末尾へ到達する全経路の確率を集約する。両者から、位置 (i,k) が対応付けられる事後確率を得る。確率の加算は非常に小さい数の underflow を避けるため log probability 上の log-sum-exp で行う。標準 beam 幅は (B_p^a=100) である。後述の復号では確率を総和せず、最も高得点の解を一つ選ぶ max-sum 計算を使う。

得られた (p^x,p^y,p^z) は `SparseFloatMatrix`、すなわち保存対象の座標と値だけを行ごとに持つ疎行列に格納し、閾値未満を保持しない。塩基対 marginal について各位置が一つの相手と対を作る確率の和は高々 1 であるため、閾値 (\theta_s) 以上の相手は一位置あたり高々 (1/\theta_s) 個である。同様に alignment match marginal の各行・列和も高々 1 なので、(\theta_a) 以上の match は一位置あたり高々 (1/\theta_a) 個である。よって固定閾値では

\[
|E_x|=O(L_1/\theta_s),\quad |E_y|=O(L_2/\theta_s),\quad
|E_z|=O((L_1+L_2)/\theta_a)
\tag{10}
\]

である。実装は sparse row の順序付き走査を用い、linear mode では decoder 入力のために密行列へ戻さない。非線形 decoder を選んだ場合に限り、互換性のため密行列を構成する。

LinearPartition または LinearAlign が内部失敗した場合は、空行列や別 model へ黙って fallback せず、非ゼロ終了で原因を通知する。これは速度比較に異なるアルゴリズムが混入することを防ぐ。

### 3.3 Profile folding と RNAalifold

profile の塩基対確率は、各構成配列で求めた BPP を現在の profile columns へ射影して平均する。以前検討した多数決 consensus sequence は compensatory substitution を失い、実測でも遅く精度が低かったため採用しない。ViennaRNA の RNAalifold は明示的な `--alifold` の場合のみ使用し、既定では無効とした。さらに folding probability engine が LinearPartition の場合は、RNAalifold が全体の線形性を破るため、指定されても無効化して配列別 BPP 平均を用いる。非線形 mode でも本研究の標準比較は `--no-alifold` であり、RNAalifold は ablation としてのみ扱う。

### 3.4 線形 max-sum folding decoder

式 (6),(7) の局所 pair score は

\[
s^x_{ij}=c^x_{ij}-q^x_{ij},\qquad
s^y_{kl}=c^y_{kl}-q^y_{kl}
\tag{11}
\]

である。`LinearNussinov` は right endpoint ごとに interval start と部分 score を持つ状態を生成し、固定幅 (B_d^s) の beam を保持する。未対合遷移、内部構造への pair 追加、保持済み prefix との結合を max で更新する。後発 start の状態が同等以上の score を持つ場合、先発 start は余分な未使用 prefix を含むだけなので suffix dominance により厳密に除ける。この dominance 除去は beam pruning ではなく、失われた導出を生じない。

確率行列に保存された候補座標と、非零ラグランジュ乗数を持つ候補座標の和集合を、塩基対の右端位置ごとに cache する。ここで cache は候補一覧を一度作って再利用する仕組みである。乗数側の候補一覧は新しい非零座標が生じたときだけ追加し、各反復で疎行列全体を再走査して並べ替え・重複除去する処理を行わない。固定 beam と閾値で疎化した候補集合の下で decoder は長さに線形である。標準 (B_d^s) は `--fold-dd-beam`、未指定時は LinearPartition beam 100 と同じである。

最終 multiple alignment の consensus structure も、平均 BPP の正の score を `LinearNussinov` で max 復号する。標準 beam (B_f^s) は `--fold-final-beam`、未指定時 100 である。

### 3.5 線形 max-sum alignment decoder

式 (8) の match score は

\[
s^z_{ik}=p^z_{ik}-\theta_a+q^z_{ik}
\tag{12}
\]

である。`LinearNeedlemanWunsch` は BeamAlign の状態遷移を probability の log-sum-exp から max へ変更し、式 (12) を match score として単調アラインメントを復号する。固定 beam (B_d^a) により線形時間・線形メモリとなる。標準値は `--align-dd-beam`、未指定時は LinearAlign probability beam 100 と同じである。

### 3.6 疎なラグランジュ乗数と一貫した subgradient

linear mode では (q^x,q^y,q^z) を疎行列で保持する。比較・デバッグ用に `--dense-lagrangian` を指定すれば密格納へ切り替えられる。CBP projection (c_x,c_y,c_z) は、各 pair または match coordinate に接続する CBP の索引であり、疎な sorted vector として保持する。

選択 CBP の個数を (n^x_{ij}=\sum_{kl}w_{ijkl})、(n^y_{kl}=\sum_{ij}w_{ijkl})、endpoint 使用数を (n^z_{ik}) とする。subgradient は、緩和した各整合性制約が現在解でどれだけ破られているかを並べたベクトルであり、乗数をどちらへ動かせば双対目的値を小さくできるかを示す。本実装では

\[
g^x_{ij}=n^x_{ij}-x_{ij},\quad
g^y_{kl}=n^y_{kl}-y_{kl},\quad
g^z_{ik}=z_{ik}-n^z_{ik}.
\tag{13}
\]

実装は非零になり得る座標を一度 `GradientEntry` という明示的な一覧に展開し、違反数、(|g|_2^2)、および実更新のすべてに同じ一覧を用いる。これにより、コンパイル時に選ばれた疎更新方式によって、更新幅の分母へ含めた座標と実際に更新した座標が食い違う問題を解消した。

### 3.7 Polyak step と二つの双対トラック

反復 (t) の実行可能下界を (LB_t)、beam 部分解が定義するラグランジアン値を (D_{\rm beam}(q_t)) とする。Polyak step は、現在の上界候補と既知の下界が離れているときは大きく、近いときは小さく乗数を動かす方法である。本実装では

\[
\eta_t=\alpha\frac{\max\{0,D_{\rm beam}(q_t)-LB_t\}}
{\|g_t\|_2^2},\qquad \alpha=0.5
\tag{14}
\]

とし、

\[
q^x\leftarrow q^x-\eta_tg^x,\quad
q^y\leftarrow q^y-\eta_tg^y,\quad
q^z\leftarrow\max\{0,q^z-\eta_tg^z\}
\tag{15}
\]

と更新する。丸め誤差で gap が負になって更新方向が反転しないよう分子を 0 で下から clamp し、ゼロ subgradient では更新しない。各 coordinate の更新絶対値は 1 に clip する。以前の独自な単調減衰は使用しない。

重要なのは、beam 解から得た (g_t) に、後述の certified relaxation 値を組み合わせないことである。beam 解はその beam objective の subgradient であり、異なる緩和関数の値を Polyak 分子へ入れると更新の整合性がない。そのため通常トラックでは (D_{\rm beam}) を式 (14) に用い、certified 値は上界記録にのみ用いる。

一方、beam 解から得た subgradient は、後述する厳密な緩和上界そのものを直接小さくする更新方向ではない。このため、固有の乗数と更新履歴を持つ第二の最適化系列を `certified` トラックとして保持する。この系列では、構造について「各左端位置から高々一塩基対」という制約だけを残す left-endpoint capacity relaxation を、alignment について「各行から高々一 match」という row-capacity relaxation を用いる。それぞれは各左端または各行の最高得点候補を選ぶだけで厳密に解け、その選択辺から対応する subgradient を得られる。certified track にも同じ Polyak 式を適用するが、分子には必ず同じ緩和関数の値を用いる。トラック数は 2 に固定されるため線形性を変えない。

### 3.8 RIBOSUM による pair-pair match のスコア化

従来は (w_{ijkl}) に固有 score がなかった。本研究では RIBOSUM85-60 を用いる。RIBOSUM は、既知の構造 RNA alignment において、ある塩基対型が別の塩基対型へ置換される頻度を背景頻度と比較して得た塩基対置換行列である。通常の nucleotide substitution matrix が一塩基対一塩基の置換を評価するのに対し、RIBOSUM は (AU) から (GC) のような二塩基からなる塩基対型同士の保存・置換を評価する。85-60 は使用する行列系列の一つを表す名称である。本実装は Infernal の公式配布行列に含まれる対称な 16×16 score を使用する。

RIBOSUM85-60 の塩基対置換 score を (R_{ijkl})、その目的関数中の重みを (\rho\ge0) とし、目的関数を

\[
\max\left\{
\sum c^x x+\sum c^y y+\sum c^z z+
\rho\sum_{(i,j,k,l)\in{\cal C}}R_{ijkl}w_{ijkl}
\right\}
\tag{16}
\]

へ拡張する。既定値は (\rho=0.075) である。

profile column pair ((i,j)) について、各構成配列で両 column が非 gap なら塩基対型 (r\in\{AA,AC,\ldots,UU\}) を数え、profile 全配列数で割った分布 (f^{(1)}_{ij}(r)) を作る。gap または未知塩基は寄与しないが、分母は非 gap 数ではなく profile size であるため、gap の多い column pair は小さく重み付けされる。二 profile の期待置換 score は

\[
R_{ijkl}=\sum_{r=1}^{16}\sum_{s=1}^{16}
f^{(1)}_{ij}(r)f^{(2)}_{kl}(s)M^{85-60}_{rs}
\tag{17}
\]

である。分布は column-pair key で cache する。

式 (5) の (w) 係数、すなわち「現在の乗数の下でこの CBP を選ぶことによるラグランジアンの増分」である reduced cost は

\[
\bar c_{ijkl}=\rho R_{ijkl}+q^x_{ij}+q^y_{kl}
-q^z_{ik}-q^z_{jl}
\tag{18}
\]

となり、(\bar c_{ijkl}>0) のとき (w_{ijkl}=1) とする。負の RIBOSUM score は pair-pair match を禁止しない。構造確率、アラインメント確率、および乗数由来の利得が負 score を上回れば選択される。したがって RIBOSUM は hard filter ではなく soft objective である。整数計画経路でも (w) 変数に同じ係数を設定し、双対経路と目的関数を一致させた。

### 3.9 動的 CBP の候補分離と候補条件

全 CBP の事前列挙を避けるため、`--dynamic-cbp` では、その時点で明示的に保持する作業候補集合 ({\cal C}_t) を反復的に更新する。まず、folding 解 (x) の各 pair ((i,j)) の両端を現在の (z) で ((k,l)) へ写し、有効なら新しい CBP 変数として追加する。逆方向にも、(y) の pair ((k,l)) を (z^{-1}) で ((i,j)) へ写す。前向きだけでは (y) にのみ現れる構造不一致を見逃すため、両方向から候補を見つける必要がある。この現在解に基づく候補発見を候補分離と呼ぶ。

候補は四つの posterior が support 上にあり、従来互換の (\rho=0) では

\[
p-\theta_s>0,\qquad
\omega(p-\theta_s)+(q-\theta_a)>0,
\tag{19}
\]

を満たすものとする。ここで

\[
p=\frac{N_1p^x_{ij}+N_2p^y_{kl}}{N_1+N_2},\qquad
q=\frac{p^z_{ik}+p^z_{jl}}{2}.
\]

RIBOSUM 使用時は各 profile と両 endpoint の objective contribution を数えた

\[
2\{\omega(p-\theta_s)+(q-\theta_a)\}+\rho R_{ijkl}>0
\tag{20}
\]

を用いる。これは負の RIBOSUM を絶対排除する条件ではない。

### 3.10 漏れのない価格付けと動的列でも失われない上界性

単に ({\cal C}_t\subset{\cal C}) だけで式 (9) を計算すると、まだ作業候補集合へ追加していない変数の正の (max(0,\bar c)) を足していないため、値は元問題の双対上界ではない。これが「動的 CBP では双対値が元問題の上界でない」という重大な問題である。

本実装は双対値を評価する前に、通常トラックまたは certified トラックのいずれかで式 (18) が正となる全ての未追加 CBP を明示的な変数として追加する。この操作が exact pricing である。ここで intrinsic score は、ラグランジュ乗数とは無関係に目的関数が CBP 自体へ与える固有得点、すなわち (\rho R_{ijkl}) を指す。正の固有得点を持つ有効列は、乗数が全てゼロでも正 reduced cost になり得るため、反復開始前に一度全て追加し、その後は削除しない。

それ以外の未追加列では (\rho R_{ijkl}\le0) かつ (q^z\ge0) である。このとき reduced cost が正なら、正の寄与を供給できる項として必ず (q^x_{ij}>0) または (q^y_{kl}>0) が存在する。従って、過去に非零になった座標も含む正の (q^x) 座標集合と正の (q^y) 座標集合の両方から、各 endpoint に許された疎な alignment match 一覧を走査すれば、正 reduced-cost 列を漏れなく発見できる。通常・certified の二系列について別々に全走査せず、両者の候補座標の和集合を一度だけ走査する。

固定閾値では一 endpoint の alignment 候補数が高々 (1/\theta_a) なので、ある structure pair から生成される四つ組は高々 (O(1/\theta_a^2)) である。よって pricing は

\[
O\left((|Q_x|+|Q_y|)/\theta_a^2\right)=O(L)
\tag{21}
\]

となる。pricing 後は ({\cal C}\setminus{\cal C}_t) の全列で (\bar c\le0) が保証され、その寄与は 0 なので、restricted sum は full (D_w(q)) と一致する。

### 3.11 長期間不活性な CBP の除去

一度追加した変数を永久保持すると、反復とともに作業候補集合が膨張する。そこで、次をすべて満たす CBP 変数を「長期間不活性」と判定する。

1. 通常・certified のいずれの (w) 解でも選択されない。
2. 通常・certified の folding 解で対応する (x_{ij}) または (y_{kl}) が選択されない。
3. 両トラックの reduced cost が非正である。
4. 上記不活性状態が連続 8 反復続く。
5. intrinsic RIBOSUM score が正でない。

alignment endpoint (z_{ik},z_{jl}) だけの選択は、consensus pair を活性化したとは数えない。条件を満たす列を除去した後、CBP の重複を高速に検出する hash set と、各構造 pair・alignment match から関連 CBP を引く索引 (c_x,c_y,c_z) を、残存列から再構築する。将来乗数が変化して reduced cost が正になれば exact pricing が再追加するため、削除は上界性を損なわない。連続 8 反復という猶予期間は、乗数が振動したときに同じ列を短期間で削除・再追加することを抑えるための実装上の定数である。

### 3.12 実行可能解への修復と単調下界

beam search が別々に返した (x,y,z,w) は、それぞれの部分問題の制約を満たしていても、相互の整合性制約を満たすとは限らない。そこで各反復で、部分問題解から元問題の全制約を満たす解を作り直す。この操作を repair と呼ぶ。(z) は decoder が返す単調 alignment なのでそのまま採用し、まず (x) と (y) が両方選び、かつ両 endpoint が (z) で対応する pair の共通部分だけを残す。pair の削除は非交差・次数制約を壊さないため、この intersection repair の解は実行可能である。

linear structure mode では、より強い consensus repair も行う。これは、二つの beam 構造の共通部分だけを使うのではなく、固定した alignment (z) の下で共通塩基対をもう一度最適化する方法である。(a) 上の候補 pair ((i,j)) を ((k,l)=(z(i),z(j))) へ写し、

\[
h_{ij}=c^x_{ij}+c^y_{kl}+\rho R_{ijkl}
\tag{22}
\]

が正で有効なものを sparse edge とする。これを `LinearNussinov` で非交差 matching として復号し、選択 pair を (b) へ写す。単調 (z) は順序を保存するため、(a) 上で非交差な集合は (b) 上でも非交差である。intersection repair と consensus repair の score の大きい方を採用する。repair decoder 自体は beam 近似でも、出力された構造は実行可能なので score は常に正しい下界である。

全反復で得た最良実行可能 score を

\[
LB_t=\max\{LB_{t-1},\operatorname{score}(\hat x_t,\hat y_t,\hat z_t,\hat w_t)\}
\tag{23}
\]

として保持し、対応する解も保存する。従って LB は単調非減少で、最終出力は最後の beam 解ではなく best feasible 解である。複数構造 threshold を用いる場合、repair は最大 threshold を用いて保守的に score を計算し、下界性を保つ。

### 3.13 folding beam に対する厳密な上界証明値

beam score は、保持された状態から見つかった構造の得点であり、(D_x,D_y) の真の最大値以下である。このため、そのまま双対上界には使えない。`LinearNussinov::decode_certified` は、beam 幅の制限によって除去した各状態を記録し、その状態を完全な構造へ伸ばしたときに得られる最大 score を別の不等式で上から抑える。この付加情報を上界証明値と呼ぶ。

各配列位置をグラフの頂点、候補塩基対を二頂点を結ぶ辺とみなし、正の pair edge score を (s_{ij}^{+}=\max(0,s_{ij})) とする。各頂点 (i) に非負の上界予算に相当する potential (u_i\ge0) を割り当て、全 edge で

\[
u_i+u_j\ge s^+_{ij}
\tag{24}
\]

を満たせば、どの塩基も高々一つの辺にしか使えないことから、任意の塩基対 matching の score は (sum_i u_i) 以下である。この条件を満たす potential 集合を feasible vertex cover と呼ぶ。実装は、各位置を左端として持つ辺の最大値、右端として持つ辺の最大値、接続辺最大値の半分、という三つの初期割当てを作る。さらに、配列を前向き・後ろ向きに計 4 回走査し、全不等式を保ったまま各 potential を小さくする。更新は

\[
u_i\leftarrow\max_{j:(i,j)\in E}(s^+_{ij}-u_j)
\tag{25}
\]

である。各割当てについて、配列先頭から各位置までの potential 累積和を保存する。これにより、枝刈りされた区間状態の既得 score に、まだ選ばれ得る位置の potential を加えた最大完成得点を定数時間で計算できる。全枝刈り状態の最大完成得点と、beam に残った全区間解の得点の最大を (U_{\rm prune})、全頂点 potential 和による上界を (U_{\rm add}) とし、

\[
U_x=\min(U_{\rm prune},U_{\rm add})
\tag{26}
\]

を folding 部分問題の上界とする。suffix dominance、すなわち「より後から始まり同等以上の score を持つ状態が存在する」という包含関係で除去された状態は、より良い保持状態で常に代替できるため証明値の対象外でよい。beam 幅による枝刈りが一度も起きなければ、beam score 自体が厳密最大である。上界は倍精度浮動小数点数で加算し、単精度へ変換するときは `nextafter` を用いて正の無限大方向へ丸める。これにより数値表現の丸めで上界を僅かに下回ることを防ぐ。

### 3.14 alignment の厳密緩和上界

単調 alignment は、二つの配列位置を左右の頂点集合とする二部グラフ上の matching、すなわち各位置を高々一回だけ対応付ける辺集合の部分集合である。単調性と第二配列側の容量制約を取り除くと、「第一配列の各位置から高々一つを選ぶ」row relaxation が得られ、

\[
U_z^{\rm row}=\sum_i\max\{0,\max_k s^z_{ik}\}.
\tag{27}
\]

同様に column relaxation は

\[
U_z^{\rm col}=\sum_k\max\{0,\max_i s^z_{ik}\}.
\tag{28}
\]

であり、

\[
U_z=\min(U_z^{\rm row},U_z^{\rm col})
\tag{29}
\]

は厳密上界である。いずれも保存された match 候補を一回走査し、各行または各列の最大正得点を足すだけなので線形時間で求まる。certified dual track では乗数更新方向も必要なため row relaxation を厳密に解き、各行で最大となった正得点 edge を返す。

### 3.15 二種類の厳密双対上界

通常トラックでは、beam 値のうち folding を (U_x,U_y)、alignment を (U_z) で置き換え、exact pricing 後の (D_w) を加える：

\[
U_{\rm cert}(q)=U_x(q)+U_y(q)+U_z(q)+D_w(q).
\tag{30}
\]

第二トラックでは、構造の left-endpoint relaxation

\[
U_x^{\rm left}=\sum_i\max\{0,\max_j s^x_{ij}\}
\tag{31}
\]

と (y) の同型式、alignment row relaxation、および同じ exact (D_w) を用いる。この関数は convex な support function の和であり、その optimizer から得た一貫した subgradient で直接最小化する。反復を通じた best upper bound は

\[
UB_t=\min\{UB_{t-1},U_{\rm cert}(q_t),
U_{\rm track}(\tilde q_t)\}
\tag{32}
\]

で、単調非増加である。停止 gap は (UB_t-LB_t)、相対 gap の報告には score scale に応じた正規化を用いる。絶対 gap が (10^{-4}\max(1,|LB_t|)) 以下なら certified gap 収束とする。

beam 解の違反がゼロでも、beam 外により良い部分解があり得るため元問題の最適性は証明されない。この場合は `agreement` ではなく `beam_stationary` として停止理由を区別する。非線形の厳密 decoder で全整合した場合のみ agreement を最適性の根拠にできる。

### 3.16 アルゴリズム全体

**Algorithm 1: 一つの profile merge に対する linear certified DAFS**

1. LinearPartition で各配列の BPP を計算し、profile columns へ射影・平均して (p^x,p^y) を得る。
2. LinearAlign で profile 間の (p^z) を得る。
3. (\theta_s,\theta_a) で疎化し、posterior support、alignment envelope、decoder support cache を作る。
4. (q=0,\tilde q=0)、(LB=-\infty,UB=+\infty)、空の動的 CBP working set を初期化する。
5. (\rho R_{ijkl}>0) の有効 CBP を一度追加する。
6. (t=0,\ldots,T-1) について：
   1. 通常トラックの (x,y,z) を固定 beam max-sum decoder で解き、folding certificate と alignment relaxation を計算する。
   2. certified track の left/row relaxation を厳密に解く。
   3. 通常解と certified 解から双方向 heuristic CBP を追加する。
   4. 両トラックの正 reduced-cost 欠落列を exact pricing で追加する。
   5. 各 CBP の式 (18) を評価し、正の (w) を選び、二つの厳密双対値を計算する。
   6. intersection/consensus repair を行い、best (LB) と解を更新する。
   7. 式 (13)–(15) で通常・certified の各乗数を更新する。
   8. 8 反復 stale な非正 intrinsic CBP を除去し、projection を再構築する。
   9. best (UB) を更新し、certified gap または stationary 条件を判定する。
7. best feasible (x,y,z) を返し、guide tree の次の profile merge へ進む。

## 4. 計算量

profile merge 一回について、(E_s=|E_x|+|E_y|)、(E_a=|E_z|)、working-set CBP 数を (C) とする。確率閾値と beam 幅を明示すると概略は次の通りである。

| 処理 | 時間 | メモリ | 備考 |
|---|---:|---:|---|
| LinearPartition | (O(B_p^s L))（実装モデルの beam 定数を含む） | (O(B_p^s L)) | 近似 partition/BPP |
| LinearAlign posterior | (O(B_p^a L))（beam 定数を含む） | (O(B_p^a L)) | forward/backward |
| sparse posterior | (O(E_s+E_a)) | (O(E_s+E_a)) | (E_s,E_a=O(L)) |
| LinearNussinov ×2 | (O(B_d^sE_s+B_d^{s,2}L)) 程度 | (O(B_d^sL+E_s)) | 固定 beam で (O(L)) |
| LinearNeedleman–Wunsch | (O(B_d^aL)) 程度 | (O(B_d^aL)) | 固定 beam |
| folding certificate | (O(E_s+L)) | (O(E_s+L)) | cover sweep 数 4 |
| alignment bound | (O(E_a+L)) | (O(L)) | row/column maxima |
| sparse gradient | (O(C+E_s+E_a)) | 同左 | projection 上のみ |
| CBP pricing | (O((|Q_x|+|Q_y|)/\theta_a^2)) | (O(C)) | 固定閾値で線形 |
| stale pruning | (O(C)) | (O(C)) | patience 8 |
| feasible repair | (O(E_s+B_d^{s,2}L)) | (O(E_s+B_d^sL)) | 線形 mode のみ強化 |

固定 (\theta_s,\theta_a,B,T) で support と (C) が上記 bounds に従えば、双対分解全体は (O(TL)=O(L)) 時間、(O(L)) メモリである。より明示的には CBP の初期・pricing scan の係数に (1/(\theta_s\theta_a^2)) が含まれる。閾値を長さとともにゼロへ近づける、beam 幅や (T) を長さに応じて増やす、または配列数を増やして全配列対 posterior を計算する場合、この線形性の主張はそのまま適用できない。

非線形 mode は従来 decoder と確率 model を残し、線形 mode との精度・速度比較および小規模な基準解に使用できる。疎/密ラグランジュ格納も切替可能である。

## 5. 実装と検証

実装は C++17 の DAFS コードベースに行った。LinearPartition-C/V、BeamAlign、疎行列、linear max decoders、`GradientManager`、relaxed bounds、RIBOSUM profile score を統合した。IP solver には利用可能なビルドで HiGHS を使える。`--max-iter=0` は双対反復なしという意味で曖昧に落とさず IP 経路を用い、HiGHS build での crash を回帰テストした。負の `--max-iter` は拒否する。

テストは以下を含む。

- LinearPartition/LinearAlign の失敗を確実に呼出元へ伝える fail-fast test。
- 小規模 exact Nussinov と linear certificate の比較、およびランダム sparse instance で `beam score <= exact optimum <= certified UB` を確認する test。
- cached support の有無で decoder/certificate が変わらないことの test。
- structure/alignment relaxation と外向き丸めの test。
- RIBOSUM matrix の対称性、既知 score、profile 期待値の test。
- linear pipeline、個別 beam option、`--no-alifold`、既定 (\rho=0.075)、既定 RNAalifold 無効、相反する option 拒否の smoke test。
- benchmark scorer の SPS、MCC、CBP 指標 test。

実行時には、一行を一つの JSON object とする機械可読形式 JSONL へ内部状態を記録する。記録項目は、入力長、model、閾値、beam、疎行列・動的 CBP の使用有無、各確率計算段階の時間と非零要素数、各 progressive merge と双対反復の beam 値、証明付き上界、best UB、LB、gap、制約違反数、Polyak 更新幅、非零乗数数、CBP の追加・価格付け・削除数、および二種類の repair score である。

benchmark runner は実験条件を読み込み、各入力を実行して結果を集約する補助プログラムである。各 merge で `LB <= UB`、best UB 非増加、LB 非減少という、正しい実装なら常に成立すべき条件を独立に検証する。予測結果、標準エラー出力、GNU `time` が測定した経過時間・CPU 時間・最大常駐メモリ、コマンドと入力の checksum、Git revision、実行ファイル checksum、計算機環境、job scheduler 情報を run ごとに保存する。checksum はファイル内容から得る固定長の識別値であり、同じ入力・実行ファイルを用いたか確認するために使う。中断後は完了済み run を識別して再利用できる。

## 6. 実験

### 6.1 条件

最新の大規模実験 `factorial-004` では、Murlet 由来 13 cases を各 3 反復、非ウイルス長鎖データの 16S 25 groups と 23S 5 groups を各 1 反復実行した。全条件で `OMP_NUM_THREADS=1`、乱数 seed 42、dynamic CBP、RNAalifold 無効、双対反復は各 progressive merge 20 とした。四条件は次の通りである。

- nonlinear / RIBOSUM off: CONTRAlign + CONTRAfold、(\rho=0)
- nonlinear / RIBOSUM on: CONTRAlign + CONTRAfold、(\rho=0.075)
- linear / RIBOSUM off: LinearAlign + LinearPartition-C (`lpc`)、(\rho=0)
- linear / RIBOSUM on: LinearAlign + LinearPartition-C、(\rho=0.075)

精度は sum-of-pairs score（SPS）、alignment PPV、structure sensitivity/PPV/MCC、および pair-pair CBP precision/recall/F1 で評価した。SPS は参照 alignment で同じ column にある塩基対のうち予測でも対応した割合、alignment PPV は予測した residue match のうち参照にも存在する割合である。structure sensitivity は参照塩基対の回収率、structure PPV は予測塩基対の正解率、MCC は正例・負例の両方を考慮する相関係数である。CBP precision/recall/F1 は、塩基対の両端が対応する pair-pair match を単位として計算し、F1 は precision と recall の調和平均である。

SCI（structure conservation index）は、multiple alignment における共通構造の安定性を個別配列の安定性と比較する指標であるが、ViennaRNA の version によって値が変わり得る。別途 version を固定する必要があるため現段階では `null` とし、環境依存の値を暗黙に混入させていない。長鎖 16S/23S には本ドラフト時点で完全な参照精度を集計しておらず、主に時間・メモリ・gap を評価した。

### 6.2 Murlet における精度と RIBOSUM 効果

| 構成 | RIBOSUM | SPS | MCC | CBP F1 | wall time (s) | max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| nonlinear | 0 | 0.8191 | 0.6784 | 0.5919 | 1.113 | 8.19 |
| nonlinear | 0.075 | 0.8224 | 0.6861 | 0.5974 | 1.124 | 8.11 |
| linear | 0 | 0.8220 | 0.6502 | 0.5690 | 2.637 | 12.84 |
| linear | 0.075 | **0.8286** | 0.6627 | 0.5828 | 2.666 | 12.92 |

値は 13 cases × 3 repetitions の run 平均である。RIBOSUM 0.075 は nonlinear で SPS +0.00325、MCC +0.00770、CBP F1 +0.00553、linear で SPS +0.00657、MCC +0.01259、CBP F1 +0.01381 を示した。計算時間増加は約 1%であった。linear は nonlinear より SPS が高い一方、構造 MCC と CBP F1 は低く、確率 model と beam 近似の差が主に構造精度へ残っていることを示す。

事前の nonlinear/no-alifold RIBOSUM parameter sweep、すなわち他の条件を固定して RIBOSUM 重みだけを順に変える実験では、(\rho=0,0.075,0.1) について、0 の SPS/MCC/F1 が 0.8191/0.6784/0.5919、0.075 が 0.8224/0.6861/0.5974、0.1 が 0.8124/0.6869/0.5958 であった。0.1 は MCC のみ僅かに高いが SPS と CBP F1 を落とすため、総合的な既定値を 0.075 とした。ただし、この選択と同じデータで最終効果を報告しているため、独立 test set による再検証が必要である。

### 6.3 長鎖 RNA の実行時間とメモリ

| データ | 構成 | RIBOSUM | wall time (s) | max RSS (MB) |
|---|---|---:|---:|---:|
| 16S (25 groups) | nonlinear | 0 | 74.40 | 219.57 |
| 16S | linear | 0 | 40.99 | 176.27 |
| 16S | linear | 0.075 | 41.22 | 177.59 |
| 23S (5 groups) | nonlinear | 0 | 517.02 | 828.56 |
| 23S | linear | 0 | 112.17 | 433.12 |
| 23S | linear | 0.075 | 112.06 | 435.75 |

RIBOSUM off 同士で、linear は no-alifold nonlinear より 16S で 1.82 倍、23S で 4.61 倍高速であった。メモリは 16S で約 19.7%、23S で約 47.7%少なかった。配列長が長い 23S ほど差が大きく、非線形処理の高次項を置換した効果と整合する。ただし、正式な線形性の実証には、同一生成分布または同一 family で長さだけを広く変え、wall time と RSS を (L,L^2,L^3) model へ回帰する追加実験が必要である。

RNAalifold を有効にしていた以前の nonlinear 実験との参考比較では、無効化により 16S が 180.14 s から 74.40 s、23S が 1276.21 s から 517.02 s へ短縮し、いずれも約 2.4 倍であった。このため標準比較で RNAalifold を無効にしたことは、linear mode だけを有利にする設定ではなく、両構成から共通の非線形 profile folding cost を除いた比較である。

### 6.4 primal–dual gap

最終相対 gap の merge 平均は次の通りであった。

| データ | linear, (\rho=0) | linear, (\rho=0.075) | nonlinear, (\rho=0) | nonlinear, (\rho=0.075) |
|---|---:|---:|---:|---:|
| Murlet | 5.84% | 5.49% | 6.21% | 5.47% |
| 16S | 2.38% | 3.13% | 1.60% | 2.37% |
| 23S | 3.27% | 3.98% | 3.41% | 4.02% |

Murlet と 23S では linear certificate の gap は nonlinear と同程度であり、16S の RIBOSUM off では約 0.78 percentage point 大きい。RIBOSUM は目的関数を強化する一方で、独立に選べる (w) の正 reduced-cost 項を増やし、緩和を弱くする場合がある。そのため精度改善と gap 縮小は必ずしも同方向ではない。

最新実装では、各反復で最も頻繁に実行される上界計算経路から、長鎖データで上界を実際には改善しなかった汎用構造緩和の再構築を外し、LinearNussinov が復号と同時に生成する証明値を再利用した。直前版 `factorial-003` に対し予測を変えず、linear/RIBOSUM-off の双対分解時間を Murlet 6.10%、16S 6.91%、23S 6.59%削減し、プログラム開始から終了までの経過時間をそれぞれ 2.05%、2.55%、3.80%削減した。

### 6.5 実験の完全性に関する注記

`factorial-004` は全 276 runs が成功し、制限時間超過も、上界・下界の単調性など常に成立すべき条件の違反もなかった。一方、実行条件一覧である manifest が記録した基点 commit は `c966081` で、benchmark 時の作業ディレクトリには、まだ commit されていないが後に `9360ad8` として commit された変更が含まれていた。実行ファイルの SHA-256 checksum は保存されているため使用バイナリ自体は同定できるが、最終論文の実験は未記録変更のない commit または release tag から再実行すべきである。また、一部の実験設定ファイルと dataset assets は当時 Git の追跡対象外であった。公開用成果物では実験条件一覧、取得元、license、checksum を固定する必要がある。

## 7. 考察

### 7.1 線形化と精度の trade-off

短い Murlet では beam と上界証明値の計算に必要な固定処理により linear 構成が nonlinear より遅かった。これは線形時間アルゴリズムが全ての有限長で高速という意味ではない。一方、16S、特に 23S では時間・メモリの優位性が明瞭になった。構造 MCC の低下は、LinearPartition-C の近似 posterior、posterior threshold、linear folding beam、LinearAlign posterior、alignment beam、progressive 誤差の複合結果である。原因切り分けのため、一度に一つの構成要素だけを従来版または線形版へ切り替える ablation experiment を行い、folding probability、alignment probability、双対分解 decoder のどの段が精度差を支配するかを定量化する必要がある。

改善候補は、(i) 確率計算用 beam と復号用 beam を別々に変える比較、(ii) (\theta_s,\theta_a) を変えて、一方の指標を改善するには他方を悪化させなければならない精度–速度の非劣解集合を調べること、(iii) RNA family または長さに応じた beam 自動選択、(iv) 疎化後 posterior の確率較正、(v) certified gap が大きい profile merge だけ beam または反復数を増やし、任意の時点で現在の解と bound を返せる制御である。ただし長さ依存で無制限に beam を増やす設定は理論上の固定 beam 線形性とは区別して報告する。

### 7.2 厳密上界の意味

本研究の UB は、LinearPartition/LinearAlign が生成した閾値疎 posterior と、候補条件で定義された DAFS 最適化問題に対して厳密である。これは真の partition function、無枝刈り posterior、あるいは生物学的正解に対する誤差 bound ではない。また、folding beam certificate と row/column relaxation は完全 decoder より緩いため、gap がゼロでないことは出力が悪いことを意味せず、単に現在の relaxation では最適性を証明できないことを意味する。

独立 certified dual track は「beam 由来 subgradient で別の上界を最小化する」という不整合を避ける。通常トラックは良い primal recovery、certified track は convex upper bound の改善という役割分担を持つ。今後は row と column、left と right の複数の certified tracks、または feasible fractional vertex-cover optimizer から subgradient を得るより強い track が候補である。ただし track 数を固定し、各 track の一反復を sparse support に線形に保つ必要がある。

### 7.3 RIBOSUM の解釈

RIBOSUM 項は、単独塩基対の存在確率だけでは区別しにくい compensatory substitution の妥当性を、pair-pair match そのものへ与える。これは (w) が単なる制約補助変数から生物学的 score を持つ決定変数へ変わったことを意味する。負 score を hard exclusion にしない設計は、RIBOSUM と posterior model の尺度不一致や、保存された非 canonical pair を許容する。

現在の profile 分布は profile size で正規化するため gap occupancy を自然に罰するが、sequence weighting や phylogenetic redundancy は考慮しない。Henikoff weighting、tree weighting、または effective sequence count による分布推定は今後の検討事項である。重み 0.075 は Murlet tuning に基づくため、nested cross-validation または family-disjoint validation が必要である。

### 7.4 CBP working set

exact pricing により、動的列集合でも各反復の (D_w) は full candidate problem と一致する。鍵は、正 intrinsic 列を常駐させ、その他の正 reduced-cost 列が正 (q^x) または正 (q^y) support を必要とするという符号論である。将来、(q^z) の符号制約を変える、別の (w) score を追加する、または候補条件を変える場合、この exhaustive 性の証明を再確認しなければならない。

stale pruning は厳密上界性を損なわないが、削除・再 pricing の定数時間を増やし得る。patience 8 は経験的であり、CBP peak、再追加率、DD 時間を目的とする ablation が望ましい。positive RIBOSUM 列を一切削除しないため、特殊な profile でその数の定数係数が大きくなる可能性も測定すべきである。

### 7.5 残る非線形依存と適用範囲

固定配列数に対する長さ方向の線形化は達成したが、multiple alignment 全体では guide tree 用の全配列対距離・posterior、profile の構成配列走査、出力サイズなどに配列数依存がある。IPknot の多 level/pseudoknot decoder、`--alifold`、CONTRAfold/CONTRAlign、密ラグランジュ、および `--max-iter=0` の整数計画は linear mode の主張に含めない。補助 posterior ファイルを与える mode も、その生成計算量は本手法外である。

## 8. 結論

DAFS の原 MEA/IP と双対分解の意味を保ちながら、posterior 計算、復号、双対格納、CBP 列管理、primal repair、upper-bound certificate を疎な固定 beam 処理へ置換した。特に、beam 値を誤って双対上界とみなさず、pruned-state certificate と凸緩和の独立 dual track を用いたこと、動的 CBP に exact pricing を加えて full problem の上界性を回復したことが理論上の中心である。さらに pair-pair match へ RIBOSUM85-60 を導入し、小さい時間増加で Murlet の alignment/structure/CBP 精度を改善した。長鎖 16S/23S では no-alifold nonlinear 構成に対し速度・メモリの改善を確認した。今後は clean release 上での length-scaling、family-disjoint parameter validation、他の長鎖 joint alignment/folding 法との比較、および構造精度低下の engine ablation が必要である。

## 付録 A. 実装パラメータ

| option | 意味 | 現在の既定値 |
|---|---|---:|
| `-a LinearAlign` | 線形近似 match posterior | 明示指定 |
| `-s lpc` | CONTRAfold parameter の LinearPartition | 明示指定 |
| `--align-beam` | alignment probability beam | 100 |
| `--align-dd-beam` | alignment DD max decoder beam | `--align-beam` と同じ |
| `--linfold-beam` | folding probability beam | 100 |
| `--fold-dd-beam` | folding DD max decoder beam | `--linfold-beam` と同じ |
| `--fold-final-beam` | final consensus folding beam | `--linfold-beam` と同じ |
| `--fold-th` | BPP threshold (\theta_s) | 0.2 |
| `--align-th` | match threshold (\theta_a) | 0.01 |
| `--weight` | structure weight (\omega) | 4.0 |
| `--ribosum-weight` | RIBOSUM weight (\rho) | **0.075** |
| `--eta` | Polyak scale (alpha) | 0.5 |
| `--max-iter` | DD iterations per merge | 600（benchmark は 20） |
| `--dynamic-cbp` | dynamic separation/pricing/pruning | 明示指定 |
| `--dense-lagrangian` | dense multiplier ablation | 無効 |
| `--alifold` | Vienna RNAalifold profile BPP | 無効（opt-in） |
| `--no-alifold` | RNAalifold 無効を明示 | 既定動作 |
| `--metrics-jsonl` | machine-readable trace | 無効 |

beam size に 0 以下が与えられた場合、linear complexity を保つため 1 に clamp し warning を出す。`--alifold` と `--no-alifold` の同時指定、負の RIBOSUM weight、負の max iteration は error とする。

## 付録 B. 投稿前に追加すべき解析

1. Murlet 全対象および DAFS 原論文の PKfree/PK benchmark を再現し、原 DAFS、linear DAFS、RIBOSUM ablation を family 単位で比較する。
2. 16S/23S の参照 alignment/structure が利用できる範囲で SPS、MCC、CBP F1 を計算する。
3. 長さ bins または controlled subsequences で wall/RSS を測定し、log–log slope と (O(L),O(L^2)) model の適合度を報告する。
4. (\rho) は training families で選び、family-disjoint test families で 0.075 の効果を検証する。
5. paired test と family-level bootstrap confidence interval、multiple-testing correction を追加する。
6. probability engine、DD decoder、RNAalifold、疎乗数、dynamic CBP、repair、certificate track を一つずつ切る ablation を行う。
7. beam/threshold/iteration と精度・時間・gap の Pareto curve を示す。
8. clean Git tag、container/module versions、compiler flags、TSUBAME queue/resource request、全 config/data checksum を公開する。
9. ViennaRNA/RNAfold を version pin した SCI を補助指標として追加する。
10. LinearTurboFold 等の長鎖 RNA multiple alignment/folding 法と、同一非ウイルスデータ上で比較する。

## 参考文献（ドラフト）

1. Sato K, Kato Y, Akutsu T, Asai K, Sakakibara Y. **DAFS: simultaneous aligning and folding of RNA sequences via dual decomposition.** *Bioinformatics*. 2012;28(24):3218–3224. https://doi.org/10.1093/bioinformatics/bts612
2. Zhang H, Zhang L, Mathews DH, Huang L. **LinearPartition: linear-time approximation of RNA folding partition function and base-pairing probabilities.** *Bioinformatics*. 2020. https://pmc.ncbi.nlm.nih.gov/articles/PMC7355276/
3. Li S, Zhang H, Mathews DH, Huang L. **LinearTurboFold: linear-time global prediction of conserved structures for RNA homologs with applications to SARS-CoV-2.** *PNAS*. 2021. https://pmc.ncbi.nlm.nih.gov/articles/PMC8719904/
4. Klein RJ, Eddy SR. **RSEARCH: finding homologs of single structured RNA sequences.** *BMC Bioinformatics*. 2003.（RIBOSUM の由来に関する引用候補。最終稿では使用した RIBOSUM85-60 配布物の一次出典を確認する。）
5. ViennaRNA Package および CONTRAfold/CONTRAlign/ProbCons の各一次文献は、最終稿で使用 version とともに追記する。
