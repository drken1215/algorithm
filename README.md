# 様々なアルゴリズムの実装例
データ構造や数論的アルゴリズムまで、様々な分野のアルゴリズムたちを C++17 で実装しています。  
アルゴリズム系の研究開発において計算機実験が必要になる場面や、
プログラミングコンテストに参加する場面などを想定して、
「実装例」または「ライブラリ」として使用することを念頭に置いています。


　


|分類|内容|具体例|
|---|---|---|
|**[MathNumberTheory](#mmt)**|整数論的アルゴリズム|素因数分解、最大公約数など|
|**[MathCombinatorics](#mc)**|組合せ論的アルゴリズム|modint、Nim など|
|**[MathAlgebra](#ma)**|代数的アルゴリズム|行列計算など|
|**[DataStructure](#ds)**|データ構造|Union-Find、セグメント木など|
|**[DataStructureOnTree](#dst)**|木上のクエリに答えるためのデータ構造|Euler ツアー、HL 分解など|
|**[GraphTheory](#gt)**|グラフアルゴリズム|強連結成分分解、木の直径など|
|**[GraphNetworkFlow](#gnf)**|ネットワークフローアルゴリズム|Ford-Fulkerson 法など|
|**[DP](#dp)**|定型的な動的計画法やその他の処理|いもす法、LIS、CHT など|
|**[Geometry](#ge)**|計算幾何|円の交点など|
|**[String](#st)**|文字列アルゴリズム|ローリングハッシュ、Suffix Array など|
|**[Others](#ot)**|その他|xorshift、サイコロなど|



　

<a name="mmt"></a>
# MathNumberTheory
整数論的アルゴリズムたちです

#### 約数, 倍数

- [約数列挙](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/divisor.cpp)
- [最大公約数 (Euclid の互除法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/GCD.cpp)
- [最小公倍数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/LCM.cpp)
- [拡張 Euclid の互除法](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/extended_GCD.cpp)

#### 素数

- [素数判定](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/is_prime.cpp)
- [素因数分解](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/prime_factorization.cpp)
- [確率的素数判定 (Miller-Rabin 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/is_prime_Miller_Rabin.cpp)
- [確率的素因数分解 (Pollard のロー法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/pollard_rho.cpp)
- [Euler のトーティエント関数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/euler_function.cpp)
- [原始根](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/primitive_root.cpp)
- N 以下の素数の個数 (O(N^2/3))

#### エラトステネスの篩

- [エラトステネスの篩](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Eratosthenes.cpp)
- [エラトステネスの区間篩](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Eratosthenes_segment.cpp)
- [高速素因数分解, 約数列挙, メビウス関数 (エラトステネスの篩風)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/fast_prime_factorization_eratosthenes.cpp)
- 線形篩
- アトキンの篩

#### 方程式

- [中国剰余定理](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/chinese_reminder_theorem.cpp)
- [中国剰余定理 (Garner 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/garner.cpp)
- [ペル方程式](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Pell_equation.cpp)
- [離散対数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/mod_log.cpp)
- 平方剰余 (Tonelli–Shanks 法)

#### 有理数

- [有理数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/rational_number.cpp)
- [Stern-Brocot 木](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Stern_Brocot.cpp)

#### その他

- [多倍長整数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/big_integer.cpp)
- [ガウス整数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/gauss_integer.cpp)
- [平衡三進法展開](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/power_of_three.cpp)
- [floor sum](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/floor_sum.cpp)



　

<a name="mc"></a>
# MathCombinatorics
組合せ論的アルゴリズムたちです

#### mod, 二項係数

- [modint](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod.cpp)
- [実行時に法が決まる modint](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_runtime.cpp)
- [累乗](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_power.cpp)
- [逆元](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_inverse.cpp)
- [二項係数 (オーソドックス, n<=10^7, r<=10^7, p<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient.cpp)
- [二項係数 (愚直計算, n<=10^9, r<=10^7, p<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_naive.cpp)
- [二項係数 (漸化式計算, n<=10^9, r<=10^9, p<=10^7)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_dp.cpp)
- [二項係数 (任意 mod, n<=10^7, r<=10^7, m<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_any_mod.cpp)
- [mod の値が大きいとき](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/big_mod.cpp)

#### 様々な数

- [重複組合せ](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/combination_with_repetition.cpp)
- カタラン数
- 分割数
- スターリング数
- ベル数
- ベルヌーイ数

#### ソート

- クイックソート
- マージソート
- ヒープソート
- コムソート
- Radix ソート
- 挿入ソート
- その他のソート達

#### マトロイド

- マトロイド上の Greedy 法
- マトロイド交差

#### その他

- [Nim](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/nim.cpp)
- [LIS and LDS](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/LIS_and_LDS.cpp)


　


<a name="ma"></a>
# MathAlgebra
行列計算など代数的計算に関するアルゴリズムです

#### 行列

- [行列 (基本演算)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix.cpp)
- [行列累乗, ランク, 連立一次方程式 (実数)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_double.cpp)
- [行列累乗, ランク, 連立一次方程式 (mod. p)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_modp.cpp)
- [行列累乗, ランク, 連立一次方程式 (binary)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_binary.cpp)
- Toeplitz 行列 (乗算, 連立方程式が O(n^2))
- 巡回行列 (乗算が O(n^2))
- コンパニオン行列
- 三重対角行列 (連立方程式が O(n))
- Black Box Linear Algebra

#### 多項式, 方程式

- 二次方程式
- 多項式 (実数係数)
- 多項式 (mod. p 係数)
- きたまさ法 (俗称)
- きたまさ法 with FFT (俗称)
- 多項式補間

#### 畳み込み計算

- [FFT (高速フーリエ変換)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/FFT.cpp)
- [NTT (高速剰余変換)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/NTT.cpp)
- [高速アダマール変換](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/fast_hadmard_transform.cpp)
- 高速ゼータ変換
- 高速メビウス変換
- [添字 GCD 畳み込み](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/fast_gcd_convolution.cpp)
- Karatsuba 法

#### 形式的冪級数

- [形式的冪級数](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/formal_power_series.cpp)
- [形式的冪級数 (実行時 mod)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/formal_power_series_runtime_mod.cpp)

#### 最適化

- 二分探索法 (方程式の解を 1 つ求める)
- 三分探索法
- 黄金探索法
- Newton 法
- 単体法
- 分枝限定法


　


<a name="ds"></a>
# DataStructure
各種データ構造の実装です

#### Union-Find

- [Union-Find (union by size)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_size.cpp)
- [Union-Find (union by rank)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_rank.cpp)
- [重みつき Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/weighted_union_find.cpp)
- [重みつき Union-Find (F2 体)](https://github.com/drken1215/algorithm/blob/master/DataStructure/weighted_union_find_F2.cpp)
- [部分永続 Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/partially_persistent_union_find.cpp)
- [undo つき Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/undo_union_find.cpp)
- Quick Find
- Dynamic Connectivity

#### セグメント木

- [セグメント木](https://github.com/drken1215/algorithm/blob/master/DataStructure/segment_tree.cpp)
- [セグメント木 (遅延評価)](https://github.com/drken1215/algorithm/blob/master/DataStructure/segment_tree_lazy.cpp)
- [Starry Sky 木 (俗称)](https://github.com/drken1215/algorithm/blob/master/DataStructure/starry_sky_tree.cpp)
- マージソート過程保存木
- 等差数列区間加算木
- 二次元セグメント木

#### Binary Indexed 木

- [BIT](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree.cpp)
- [BIT 上二分探索 (k 番目の要素を求める)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_search_on_BIT.cpp)
- [BIT (区間加算, 区間和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_RAQ.cpp)
- [二次元 BIT](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_2D.cpp)
- [二次元 BIT (領域加算, 領域和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_2D_RAQ.cpp)

#### RMQ

- [RMQ (セグメント木)](https://github.com/drken1215/algorithm/blob/master/DataStructure/range_minimum_query.cpp)
- [RMQ (Sparse Table)](https://github.com/drken1215/algorithm/blob/master/DataStructure/sparse_table.cpp)

#### 平方分割

- [Mo 法](https://github.com/drken1215/algorithm/blob/master/DataStructure/mo.cpp)

#### 平衡二分探索木

- [RBST](https://github.com/drken1215/algorithm/blob/master/DataStructure/randomized_binary_search_tree.cpp)
- Treap 木
- AVL 木
- Splay 木
- 赤黒木

#### 永続データ構造

- 永続配列
- 完全永続 Union-Find 木
- 永続セグメント木
- 永続赤黒木

#### ハッシュ

- [Zobrist hash](https://github.com/drken1215/algorithm/blob/master/DataStructure/zobrist_hash.cpp)
- 木に対する hash

#### ヒープ

- Skew Heap (マージ可能)
- Paring Heap (マージ可能)
- Radix Heap
- Fibonacci Heap

#### その他

- [Binary Trie](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_trie.cpp)
- [Disjoint Sparse Table](https://github.com/drken1215/algorithm/blob/master/DataStructure/disjoint_sparse_table.cpp)
- [並列二分探索](https://github.com/drken1215/algorithm/blob/master/DataStructure/parallel_binary_search.cpp)
- [Cartesian 木](https://github.com/drken1215/algorithm/blob/master/DataStructure/cartesian_tree.cpp)
- Wavelet 木

　

<a name="dst"></a>

# DataStructureOnTree
ツリー上のクエリ処理のためのデータ構造たちの実装です

#### 木全般

- [部分木サイズ, LCA など](https://github.com/drken1215/algorithm/blob/master/DataStructureOnTree/run_tree.cpp)

#### LCA

- [LCA (ダブリング)](https://github.com/drken1215/algorithm/blob/master/DataStructureOnTree/LCA_doubling.cpp)
- LCA (Euler Tour)
- LCA (HL 分解)

#### テクニック

- [Euler Tour (頂点上)](https://github.com/drken1215/algorithm/blob/master/DataStructureOnTree/euler_tour_on_nodes.cpp)
- [Euler Tour (辺上のクエリ)](https://github.com/drken1215/algorithm/blob/master/DataStructureOnTree/euler_tour_on_edges.cpp)
- [HL 分解](https://github.com/drken1215/algorithm/blob/master/DataStructureOnTree/heavy_light_decomposition.cpp)
- 重心分解
- Link-Cut 木
- マージテク (俗称)
- DSU on Tree

#### その他の問題

- Level Ancester


　

<a name="gt"></a>

# GraphTheory
グラフ理論全般のアルゴリズムです

#### DFS・BFS

- DFS (連結成分を数える)
- BFS (重みなしグラフの最短路)
- トポロジカルソート (DFS)
- トポロジカルソート (BFS)
- サイクル検出 (DFS)
- サイクル検出 (BFS)
- サイクル検出 (Union-Find)
- [二部グラフ判定 (DFS)](https://github.com/drken1215/algorithm/blob/master/GraphTheory/is_bipartite_dfs.cpp)
- 二部グラフ判定 (BFS)
- 二部グラフ判定 (Union-Find)

#### 連結成分分解

- [強連結成分分解](https://github.com/drken1215/algorithm/blob/master/GraphTheory/strongly_connected_components.cpp)
- [橋, 関節点列挙 (Low-Link)](https://github.com/drken1215/algorithm/blob/master/GraphTheory/low_link.cpp)
- [二重辺連結成分分解](https://github.com/drken1215/algorithm/blob/master/GraphTheory/two_edge_connected_components.cpp)
- 二重頂点連結成分分解
- 2-SAT

#### ツリー

- [ツリーの直径](https://github.com/drken1215/algorithm/blob/master/GraphTheory/diameter.cpp)
- ツリーの重心

#### 最短路

- 重みなしグラフの最短路 (BFS)
- 重みが 0, 1 のみのグラフの最短路 (0-1 BFS)
- 単一始点最短路 (Dijkstra 法, 正辺のみ)
- 単一始点最短路 (Bellman-Ford 法, 負辺対応)
- [全頂点対間最短路 (Floyd-Warshall 法)](https://github.com/drken1215/algorithm/blob/master/GraphTheory/floyd_warshall.cpp)
- 全頂点対間最短路 (Johnson 法)
- k-最短路
- SPFA

#### その他

- 最小全域木 (Kruskal 法)
- 最小有向全域木 (Chu-Liu/Edmonds 法)
- 有向 Euler 路
- [無向 Euler 路](https://github.com/drken1215/algorithm/blob/master/GraphTheory/euler_tour_undirected.cpp)
- [彩色数 (O(n2^n))](https://github.com/drken1215/algorithm/blob/master/GraphTheory/vertex_coloring.cpp)
- [最大安定集合問題 (O(1.381^n))](https://github.com/drken1215/algorithm/blob/master/GraphTheory/maximum_stable_set.cpp)
- 最大クリーク列挙（O(1.443^n)）
- 最小シュタイナー木 (O(n 3^t + n^2 2^t + n^3))


　


<a name="gnf"></a>
# GraphNetworkFlow
グラフネットワークフロー関連のアルゴリズムです

#### 最大流

- [最大流 (Ford-Fulkerson 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/max_flow_ford_fulkerson.cpp)
- [最大流 (Dinic 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/max_flow_dinic.cpp)

#### 最小費用流

- [最小費用流 (Primal-Dual 法, 正辺のみ)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_flow_primal_dual.cpp)
- [最小費用流 (Primal-Dual 法, 負辺対応)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_flow_primal_dual_negative.cpp)
- 最小費用最大流 (Primal-Dual 法, 正辺のみ)
- 最小費用最大流 (Primal-Dual 法, 負辺対応)
- 最小費用循環流 (Cost-Scaling, 負閉路OK)

#### カット

- 最小カット (= 最大流)
- 全域最小カット（Stoer-Wanger 法）
- 全頂点対間最小カット (Nagamochi-Ibaraki 法)
- Gomory-Hu 木

#### マッチング

- [二部マッチング (Hopcroft-Karp 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/hopcroft_karp.cpp)
- 重みつき二部マッチング (Hungarian 法)
- 一般グラフの最大マッチング (Edmonds 法)
- 一般グラフの最大マッチング (行列補間)


　


<a name="dp"></a>
# DP
定型的な動的計画法やその他の処理です

#### 典型処理

- [累積和](https://github.com/drken1215/algorithm/blob/master/DP/cumulative_sum.cpp)
- [二次元累積和](https://github.com/drken1215/algorithm/blob/master/DP/cumulative_sum_2D.cpp)
- [いもす法 (俗称)](https://github.com/drken1215/algorithm/blob/master/DP/imos.cpp)
- [二次元いもす法 (俗称)](https://github.com/drken1215/algorithm/blob/master/DP/imos_2D.cpp)
- [スライド最小値](https://github.com/drken1215/algorithm/blob/master/DP/sliding_minimum.cpp)

#### 典型的 DP

- [転倒数](https://github.com/drken1215/algorithm/blob/master/DP/inversion_number.cpp)
- [LIS](https://github.com/drken1215/algorithm/blob/master/DP/longest_increasing_sequence.cpp)
- LCS
- 編集距離
- 重みつき区間スケジューリング問題
- ヒストグラム長方形面積最大化
- 最適二分探索木
- Set Cover
- k-Cover (O(n 2^n))
- k-partition (O(n^3 2^n))

#### DP パターン例

- ナップサック DP
- 区間分割型ナップサック DP
- bitDP
- 桁 DP
- 部分列 DP
- ダブリング DP
- 木 DP
- 全方位木 DP (俗称)
- 二乗の木 DP (俗称)

#### DP 高速化テクニック

- 累積和
- スライド最小値
- インライン DP (俗称)
- [Convex Hull Trick (傾き単調, クエリも単調)](https://github.com/drken1215/algorithm/blob/master/DP/convex_hull_trick_both_monotone.cpp)
- [Convex Hull Trick (傾き単調)](https://github.com/drken1215/algorithm/blob/master/DP/convex_hull_trick_slope_monotone.cpp)
- [Convex Hull Trick (単調でなくてよい)](https://github.com/drken1215/algorithm/blob/master/DP/convex_hull_trick.cpp)
- Monotone Minima
- Divide and Conquer
- Monge
- Alien DP
- 戻す DP (俗称)


　


<a name="ge"></a>
# Geometry
幾何ライブラリです

- [全部乗せ](https://github.com/drken1215/algorithm/blob/master/Geometry/All.cpp)
- [基本要素 (点, 線分, 円)](https://github.com/drken1215/algorithm/blob/master/Geometry/basic_elements.cpp)

#### 点, 線分, 三角形などの位置関係

- [点と線分の位置関係 (ccw)](https://github.com/drken1215/algorithm/blob/master/Geometry/ccw.cpp)
- [点と三角形の包含関係](https://github.com/drken1215/algorithm/blob/master/Geometry/is_contain_in_the_triangle.cpp)

#### 射影, 交差判定, 距離

- [射影](https://github.com/drken1215/algorithm/blob/master/Geometry/projection.cpp)
- [線分と線分の交差判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_intersect_two_segments.cpp)
- [線分と線分との距離](https://github.com/drken1215/algorithm/blob/master/Geometry/distance_two_segments.cpp)

#### 直線や円の交点

- [直線と直線の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_two_lines.cpp)
- [円と直線の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_line_circle.cpp)
- [円と円の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_two_circles.cpp)
- [円と線分の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_segment_circle.cpp)

#### 多角形

- [多角形の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_polygon.cpp)
- [点と多角形の包含判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_contain_in_the_polygon.cpp)
- [凸性判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_convex.cpp)
- [凸包](https://github.com/drken1215/algorithm/blob/master/Geometry/convex_hull.cpp)
- [凸多角形の切断](https://github.com/drken1215/algorithm/blob/master/Geometry/convex_cut.cpp)
- [ボロノイ図 (単純ver, O(N^3))](https://github.com/drken1215/algorithm/blob/master/Geometry/voronoi.cpp)
- 凸多角形の直径
- [円と円の共通部分の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_common_two_circles.cpp)
- [円と多角形との共通部分の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_common_circle_polygon.cpp)

#### 接線

- [接線](https://github.com/drken1215/algorithm/blob/master/Geometry/tanline.cpp)
- [共通接線 (2 円)](https://github.com/drken1215/algorithm/blob/master/Geometry/common_tanline.cpp)

#### 三次元幾何

- [三次元幾何一式](https://github.com/drken1215/algorithm/blob/master/Geometry/basic_elements_3D.cpp)

#### その他

- [最近点対](https://github.com/drken1215/algorithm/blob/master/Geometry/closest_two_points.cpp)
- 最近円対
- 線分併合
- 線分アレンジメント
- 3 点を通る円
- アポロニウスの円
- 最小包含円
- 双対変換
- kd 木


　


<a name="st"></a>
# String
文字列アルゴリズムです

#### 構文解析

- LL(1) 再帰降下パーサ

#### 文字列検索

- [ローリングハッシュ](https://github.com/drken1215/algorithm/blob/master/String/rolling_hash.cpp)
- 二次元ローリングハッシュ
- [単一パターン検索 (KMP 法)](https://github.com/drken1215/algorithm/blob/master/String/knuth_morris_pratt.cpp)
- 単一パターン検索 (Boyer-Moore 法)
- 複数パターン検索 (Aho-Corasick 法)

#### 文字列系アルゴリズム

- [Z 法](https://github.com/drken1215/algorithm/blob/master/String/z_algorithm.cpp)
- [Manacher 法](https://github.com/drken1215/algorithm/blob/master/String/manacher.cpp)

#### 文字列系データ構造

- Trie 木
- [Suffix Array](https://github.com/drken1215/algorithm/blob/master/String/suffix_array.cpp)
- Palindromic 木 (AOJ 2292)

#### その他

- [各 index 以降で各文字が最初に登場する index を求める関数](https://github.com/drken1215/algorithm/blob/master/String/next.cpp)
- split 関数
- 二次元盤面に番兵追加
- 二次元盤面を 90 度回転


　


<a name="ot"></a>
# Others
その他のアルゴリズムです

#### グリッド

- [グリッドの 4 近傍, 8 近傍](https://github.com/drken1215/algorithm/blob/master/Others/grid_neighbors.cpp)
- [ハニカムの 6 近傍](https://github.com/drken1215/algorithm/blob/master/Others/honeycomb_neighbors.cpp)

#### ビット演算テクニック

- [XorShift, ランダムシャッフル](https://github.com/drken1215/algorithm/blob/master/Others/xorshift.cpp)
- [next_combination](https://github.com/drken1215/algorithm/blob/master/Others/next_combination.cpp)
- [部分集合の部分集合](https://github.com/drken1215/algorithm/blob/master/Others/subset_enumeration.cpp)

#### 探索法

- α-β 探索
- 焼き鈍し法
- A*
- IDA*
- Baby-Step Giant-Step 法
- 平面走査法

#### その他

- [デバッグストリーム, chmin, chmax](https://github.com/drken1215/algorithm/blob/master/Others/debug.cpp)
- [pn + r (n は非負整数) で表せる整数のうち, x 以上となる最小の整数](https://github.com/drken1215/algorithm/blob/master/Others/amari_lower_bound.cpp)
- [タイマー](https://github.com/drken1215/algorithm/blob/master/Others/timer.cpp)
- [サイコロ](https://github.com/drken1215/algorithm/blob/master/Others/dice.cpp)
- [曜日](https://github.com/drken1215/algorithm/blob/master/Others/day_of_the_week.cpp)
- 四面体 (AOJ 2060)




　



# License
These codes are licensed under CC0.
[![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png "CC0")](http://creativecommons.org/publicdomain/zero/1.0/deed.ja)

