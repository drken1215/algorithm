# 様々なアルゴリズムの実装例
データ構造や数論的アルゴリズムまで、様々な分野のアルゴリズムたちを C++17 で実装しています。  
アルゴリズム系の研究開発において計算機実験が必要になる場面や、
プログラミングコンテストに参加する場面などを想定して、
「実装例」または「ライブラリ」として使用することを念頭に置いています。

　

|分類|内容|具体例|
|---|---|---|
|**[DATA STRUCTURE](#ds)**|各種データ構造|Union-Find、Sparse Table など|
|**[DATA STRUCTURE : SEGMENT](#dss)**|区間クエリに強いデータ構造|セグメント木、BIT など|
|**[DP](#dp)**|定型的な動的計画法やその他の処理|いもす法、LIS、CHT など|
|**[GEOMETRY](#ge)**|計算幾何|円の交点など|
|**[GRAPH](#gt)**|グラフアルゴリズム|強連結成分分解など|
|**[GRAPH : NETWORK FLOW](#gnf)**|ネットワークフローアルゴリズム|Ford-Fulkerson 法など|
|**[MATH : ALGEBRA](#ma)**|代数的アルゴリズム|行列計算など|
|**[MATH : COMBINATORICS](#mc)**|組合せ論的アルゴリズム|modint、Nim など|
|**[MATH : NUMBER THEORY](#mmt)**|整数論的アルゴリズム|素因数分解、最大公約数など|
|**[STRING](#st)**|文字列アルゴリズム|Suffix Array、KMP 法など|
|**[TREE](#tr)**|木上のデータ構造とアルゴリズム|Euler ツアー、木の直径など|
|**[OTHERS](#ot)**|その他|xorshift、サイコロなど|


　
<a name="ds"></a>
# データ構造 (DATA STRUCTURE)
各種データ構造の実装です

### Union-Find

- (★☆☆☆) [Union-Find (union by size)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_size.cpp)
- (★☆☆☆) [Union-Find (union by rank)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_rank.cpp)
- (★★☆☆) [重みつき Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/weighted_union_find.cpp)
- (★★★☆) [undo つき Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_can_undo.cpp)
- (★★★☆) [部分永続 Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/partially_persistent_union_find.cpp)

### キュー, ヒープ

- (★☆☆☆) [二分ヒープ](https://github.com/drken1215/algorithm/blob/master/DataStructure/heap.cpp)
- (★★☆☆) 両端 priority queue
- (★★☆☆) 削除可能 priority queue
- (★★★☆) SWAG
- (★★★★) Skew Heap (マージ可能ヒープ)
- (★★★★) Paring Heap (マージ可能ヒープ)
- (★★★★) Radix Heap
- (★★★★) Fibonacci Heap

### ハッシュ

- (★★☆☆) [ハッシュ](https://github.com/drken1215/algorithm/blob/master/DataStructure/hash.cpp)
- (★★★☆) [Zobrist ハッシュ](https://github.com/drken1215/algorithm/blob/master/DataStructure/zobrist_hash.cpp)
- (★★★☆) ローリングハッシュ
- (★★★★) 根付き木のハッシュ

### さまざまな木

- (★★★☆) [Cartesian 木](https://github.com/drken1215/algorithm/blob/master/DataStructure/cartesian_tree.cpp)
- (★★★★) [Binary Trie](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_trie.cpp)
- (★★★★) van Emde Boas 木
- (★★★★) 64 分木

### 永続データ構造

- (★★★★) 永続配列
- (★★★★) 完全永続 Union-Find
- (★★★★) 永続キュー
- (★★★★) 永続セグメント木
- (★★★★) 永続赤黒木

### 各種高速化アルゴリズム

- (★☆☆☆) [累積和](https://github.com/drken1215/algorithm/blob/master/DP/cumulative_sum.cpp)
- (★☆☆☆) [二次元累積和](https://github.com/drken1215/algorithm/blob/master/DP/cumulative_sum_2D.cpp)
- (★☆☆☆) [いもす法 (俗称)](https://github.com/drken1215/algorithm/blob/master/DP/imos.cpp)
- (★★☆☆) [二次元いもす法 (俗称)](https://github.com/drken1215/algorithm/blob/master/DP/imos_2D.cpp)
- (★★★☆) [スライド最小値](https://github.com/drken1215/algorithm/blob/master/DP/sliding_minimum.cpp)
- (★★★☆) [Mo 法](https://github.com/drken1215/algorithm/blob/master/DataStructure/mo.cpp)
- (★★★☆) [並列二分探索](https://github.com/drken1215/algorithm/blob/master/DataStructure/parallel_binary_search.cpp)

### その他

- (★★★☆) [区間の集合を set で管理する](https://github.com/drken1215/algorithm/blob/master/DataStructure/intervals_management.cpp)
- (★★★★) Dynamic Connectivity


　
<a name="dss"></a>
# 区間系データ構造 (DATA STRUCTURE : SEGMENT)
セグメント木や BIT など、区間クエリに強いデータ構造の実装です

### セグメント木

- (★★☆☆) [セグメント木](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/segment_tree.cpp)
- (★★★☆) [セグメント木 (遅延評価)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/segment_tree_lazy.cpp)
- (★★☆☆) [RMQ (セグメント木)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/range_minimum_query.cpp)
- (★★★☆) [Starry Sky 木 (俗称)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/starry_sky_tree.cpp)
- (★★★★) Segment Tree Beats (俗称)

### Binary Indexed 木

- (★★☆☆) [BIT](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_indexed_tree.cpp)
- (★★★☆) [BIT 上二分探索 (k 番目の要素を求める)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_search_on_BIT.cpp)
- (★★★☆) [BIT (区間加算, 区間和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_indexed_tree_RAQ.cpp)

### 二次元セグメント木

- (★★★☆) 二次元セグメント木
- (★★★☆) [二次元 BIT](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_indexed_tree_2D.cpp)
- (★★★★) [二次元 BIT (領域加算, 領域和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_indexed_tree_2D_RAQ.cpp)
- (★★★★) 動的二次元セグメント木
- (★★★★) 動的二次元 BIT

### Sparse Table

- (★★★☆) [Sparse Table](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/sparse_table.cpp)
- (★★★★) [Disjoint Sparse Table](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/disjoint_sparse_table.cpp)

### ウェーブレット行列

- (★★★★) [ウェーブレット行列](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/wavelet_matrix.cpp)
- (★★★★) [BIT on ウェーブレット行列 (一点加算対応)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/BIT_on_wavelet_matrix.cpp)
- (★★★★) 動的ウェーブレット行列

### 平衡二分探索木

- (★★★★) [RBST](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/randomized_binary_search_tree.cpp)
- (★★★★) Treap 木
- (★★★★) AVL 木
- (★★★★) Splay 木
- (★★★★) 赤黒木


　
<a name="dp"></a>
# 動的計画法 (DP)
定型的な動的計画法やその他の処理です

### 有名問題

- (★☆☆☆) [ナップサック問題](https://github.com/drken1215/algorithm/blob/master/DP/knapsack.cpp)
- (★☆☆☆) [LIS](https://github.com/drken1215/algorithm/blob/master/DP/longest_increasing_sequence.cpp)
- (★☆☆☆) [LCS](https://github.com/drken1215/algorithm/blob/master/DP/lcs.cpp)
- (★★☆☆) TSP (O(2^N N^2))
- (★★☆☆) 部分列の個数
- (★★★☆) 最適二分探索木
- (★★★☆) Set Cover
- (★★★☆) k-Cover (O(2^N N))
- (★★★☆) k-partition (O(2^N N^3))

### グリッドに関する DP

- (★☆☆☆) [編集距離](https://github.com/drken1215/algorithm/blob/master/DP/edit_distance.cpp)
- (★★☆☆) [グリッドに含まれる最大正方形](https://github.com/drken1215/algorithm/blob/master/DP/largest_square_in_grid.cpp)
- (★★★☆) [ヒストグラム長方形面積最大化](https://github.com/drken1215/algorithm/blob/master/DP/histogram.cpp)

### Convex Hull Trick

- (★★★☆) [Convex Hull Trick (傾き単調, クエリも単調)](https://github.com/drken1215/algorithm/blob/master/DP/convex_hull_trick_both_monotone.cpp)
- (★★★☆) [Convex Hull Trick (傾き単調)](https://github.com/drken1215/algorithm/blob/master/DP/convex_hull_trick_slope_monotone.cpp)
- (★★★★) [Convex Hull Trick (単調でなくてよい)](https://github.com/drken1215/algorithm/blob/master/DP/convex_hull_trick.cpp)

### DP 高速化技法

- (★★★★) Slope Trick
- (★★★★) Monotone Minima
- (★★★★) Alien DP
- (★★★★) LARSCH


　
<a name="ge"></a>
# 幾何 (GEOMETRY)
幾何ライブラリです

### 全般

- (★★★☆) [全部乗せ](https://github.com/drken1215/algorithm/blob/master/Geometry/All.cpp)
- (★★☆☆) [基本要素 (点, 線分, 円)](https://github.com/drken1215/algorithm/blob/master/Geometry/basic_elements.cpp)
- (★★☆☆) [偏角ソート](https://github.com/drken1215/algorithm/blob/master/Geometry/arg_sort.cpp)

### 点, 線分, 三角形などの位置関係

- (★★☆☆) [点と線分の位置関係 (ccw)](https://github.com/drken1215/algorithm/blob/master/Geometry/ccw.cpp)
- (★★☆☆) [点と三角形の包含関係](https://github.com/drken1215/algorithm/blob/master/Geometry/is_contain_in_the_triangle.cpp)

### 射影, 交差判定, 距離

- (★★☆☆) [射影](https://github.com/drken1215/algorithm/blob/master/Geometry/projection.cpp)
- (★★☆☆) [線分と線分の交差判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_intersect_two_segments.cpp)
- (★★☆☆) [線分と線分との距離](https://github.com/drken1215/algorithm/blob/master/Geometry/distance_two_segments.cpp)

### 直線や円の交点

- (★★☆☆) [直線と直線の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_two_lines.cpp)
- (★★★☆) [円と直線の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_line_circle.cpp)
- (★★★☆) [円と円の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_two_circles.cpp)
- (★★★☆) [円と線分の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_segment_circle.cpp)

### 接線

- (★★★☆) [接線](https://github.com/drken1215/algorithm/blob/master/Geometry/tanline.cpp)
- (★★★☆) [共通接線 (2 円)](https://github.com/drken1215/algorithm/blob/master/Geometry/common_tanline.cpp)

### 多角形

- (★★☆☆) [多角形の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_polygon.cpp)
- (★★☆☆) [点と多角形の包含判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_contain_in_the_polygon.cpp)
- (★★☆☆) [凸性判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_convex.cpp)
- (★★★☆) [凸包](https://github.com/drken1215/algorithm/blob/master/Geometry/convex_hull.cpp)
- (★★★☆) [凸多角形の直径](https://github.com/drken1215/algorithm/blob/master/Geometry/diameter.cpp)
- (★★★☆) [凸多角形の切断](https://github.com/drken1215/algorithm/blob/master/Geometry/convex_cut.cpp)
- (★★★☆) [ボロノイ図 (単純ver, O(N^3))](https://github.com/drken1215/algorithm/blob/master/Geometry/voronoi.cpp)
- (★★★★) [円と円の共通部分の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_common_two_circles.cpp)
- (★★★★) [円と多角形との共通部分の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_common_circle_polygon.cpp)

### 三次元幾何

- (★★★★) [三次元幾何一式](https://github.com/drken1215/algorithm/blob/master/Geometry/basic_elements_3D.cpp)

### その他

- (★★★☆) [最近点対](https://github.com/drken1215/algorithm/blob/master/Geometry/closest_two_points.cpp)
- (★★★☆) 線分併合
- (★★★☆) 線分アレンジメント
- (★★★☆) 3 点を通る円
- (★★★☆) [アポロニウスの円](https://github.com/drken1215/algorithm/blob/master/Geometry/apollonius.cpp)
- (★★★☆) 双対変換
- (★★★★) 最小包含円
- (★★★★) kd 木


　
<a name="gt"></a>
# グラフ (GRAPH)
グラフアルゴリズムです

### DFS, BFS

- (★☆☆☆) [連結成分の個数 (DFS)](https://github.com/drken1215/algorithm/blob/master/Graph/dfs.cpp)
- (★☆☆☆) [連結成分の個数 (BFS)](https://github.com/drken1215/algorithm/blob/master/Graph/bfs.cpp)
- (★☆☆☆) [二部グラフ判定 (DFS)](https://github.com/drken1215/algorithm/blob/master/Graph/is_bipartite_dfs.cpp)
- (★☆☆☆) [二部グラフ判定 (BFS)](https://github.com/drken1215/algorithm/blob/master/Graph/is_bipartite_bfs.cpp)
- (★★☆☆) [トポロジカルソート (DFS)](https://github.com/drken1215/algorithm/blob/master/Graph/topological_sort_dfs.cpp)
- (★★☆☆) [トポロジカルソート (BFS)](https://github.com/drken1215/algorithm/blob/master/Graph/topological_sort_bfs.cpp)

### 連結成分分解

- (★★☆☆) [閉路検出 (サイクル検出)](https://github.com/drken1215/algorithm/blob/master/Graph/cycle_detection.cpp)
- (★★☆☆) [Functional グラフの閉路検出 (サイクル検出)](https://github.com/drken1215/algorithm/blob/master/Graph/functional_graph_cycle_detection.cpp)
- (★★★☆) [強連結成分分解](https://github.com/drken1215/algorithm/blob/master/Graph/strongly_connected_components.cpp)
- (★★★☆) [橋, 関節点列挙 (Low-Link)](https://github.com/drken1215/algorithm/blob/master/Graph/low_link.cpp)
- (★★★★) [二重辺連結成分分解 (Bridge-Block 木)](https://github.com/drken1215/algorithm/blob/master/Graph/two_edge_connected_components.cpp)
- (★★★★) 二重頂点連結成分分解 (Block-Cut 木)
- (★★★★) 三重辺連結成分分解 (SPQR 木)

### 最短路問題 (基本)

- (★☆☆☆) [重みなしグラフの最短路 (BFS, O(E))](https://github.com/drken1215/algorithm/blob/master/Graph/shortest_path_bfs.cpp)
- (★★☆☆) [重みが 0, 1 のみのグラフの最短路 (0-1 BFS, O(E))](https://github.com/drken1215/algorithm/blob/master/Graph/shortest_path_01bfs.cpp)
- (★★☆☆) [単一始点最短路 (Dijkstra 法, 正辺のみ, O(V + E log V))](https://github.com/drken1215/algorithm/blob/master/Graph/shortest_path_dijkstra.cpp)
- (★★☆☆) [単一始点最短路 (Bellman-Ford 法, 負辺対応, O(VE))](https://github.com/drken1215/algorithm/blob/master/Graph/shortest_path_bellman_ford.cpp)
- (★★☆☆) [全頂点対間最短路 (Floyd-Warshall 法, O(V^3))](https://github.com/drken1215/algorithm/blob/master/Graph/floyd_warshall.cpp)
- (★★★☆) 全頂点対間最短路 (Johnson 法, O(EV log V))

### 最短路問題 (応用)

- (★★★☆) SPFA
- (★★★★) 補グラフの最短路
- (★★★★) d-辺最短路
- (★★★★) Monge グラフ上の d-辺最短路

### 全域木, 路に関する問題

- (★★☆☆) 最小全域木 (Kruskal 法)
- (★★★☆) 有向 Euler 路
- (★★★☆) [無向 Euler 路](https://github.com/drken1215/algorithm/blob/master/Graph/euler_tour_undirected.cpp)

- (★★★★) 最小有向全域木 (Chu-Liu/Edmonds 法)
- (★★★★) 最小シュタイナー木 (O(V 3^t + V^2 2^t + V^3))

### グラフ上の有名問題

- (★★★☆) [最大安定集合問題 (O(1.381^V))](https://github.com/drken1215/algorithm/blob/master/Graph/maximum_stable_set.cpp)
- (★★★★) 最大クリーク列挙（O(1.443^V)）
- (★★★★) [頂点彩色 (O(2^V V))](https://github.com/drken1215/algorithm/blob/master/Graph/vertex_coloring.cpp)
- (★★★★) 辺彩色
- (★★★★) 二部グラフの辺彩色 (Alon 法, O(E log E))


　
<a name="gnf"></a>
# ネットワークフロー (GRAPH : NETWORK FLOW)
グラフネットワークフロー関連のアルゴリズムです

### 最大流

- (★★★☆) [最大流 (Ford-Fulkerson 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/max_flow_ford_fulkerson.cpp)
- (★★★☆) [最大流 (Dinic 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/max_flow_dinic.cpp)

### 最小費用流

- (★★★☆) [最小費用流 (Primal-Dual 法, 正辺のみ)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_flow_primal_dual.cpp)
- (★★★☆) [最小費用流 (Primal-Dual 法, 負辺対応)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_flow_primal_dual_negative.cpp)
- (★★★★) [最小費用循環流 (Goleberg-Tarjan 法, by cost-scaling, 負閉路 OK)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_circulating_flow.cpp)

### カット

- (★★★☆) 最小カット (= 最大流)
- (★★★★) 全域最小カット（Stoer-Wanger 法）
- (★★★★) 全頂点対間最小カット (Nagamochi-Ibaraki 法)
- (★★★★) Gomory-Hu 木

### マッチング

- (★★★☆) [二部マッチング (Hopcroft-Karp 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/hopcroft_karp.cpp)
- (★★★☆) 重みつき二部マッチング (Hungarian 法)
- (★★★★) 一般グラフの最大マッチング (Edmonds 法)
- (★★★★) 一般グラフの最大マッチング (行列補間)
- (★★★★) 重み付き一般グラフの最大マッチング


　
<a name="ma"></a>
# 代数 (MATH : ALGEBRA)
行列計算など代数的計算に関するアルゴリズムです

### 行列

- (★★☆☆) [行列 (基本演算)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix.cpp)
- (★★★☆) [行列累乗, ランク, 連立一次方程式 (実数)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_double.cpp)
- (★★★☆) [行列累乗, ランク, 連立一次方程式 (mod. p)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_modp.cpp)
- (★★★☆) [行列累乗, ランク, 連立一次方程式 (binary)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_binary.cpp)
- (★★★★) Black Box Linear Algebra

### 特殊な行列

- (★★★★) Toeplitz 行列 (乗算, 連立方程式が O(n^2))
- (★★★★) 巡回行列 (乗算が O(n^2))
- (★★★★) コンパニオン行列
- (★★★★) 三重対角行列 (連立方程式が O(n))

### FFT, NTT, Convolution

- (★★★☆) [FFT (高速フーリエ変換)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/FFT.cpp)
- (★★★☆) [NTT (高速剰余変換)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/NTT.cpp)
- (★★★★) [任意 mod 畳み込み](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/NTT_any_mod.cpp)
- (★★★☆) [添字 GCD 畳み込み](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/fast_gcd_convolution.cpp)
- (★★★★) Karatsuba 法

### 多項式 (Polynomial)

- (★★★★) [多項式：全部乗せ](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial.cpp)
- (★★★☆) [多項式の乗算 (by NTT, O(N log N))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_mul.cpp)
- (★★★★) [多項式の除算 (by NTT, inv of FPS, O(N log N))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_div.cpp)
- (★★★★) 多項式補間
- (★★★☆) [Polynomial Taylor Shift](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_taylor_shift.cpp)
- (★★★★) 多項式 GCD

### 形式的冪級数 (FPS)

- (★★★★) [形式的冪級数：全部乗せ](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/formal_power_series.cpp)
- (★★★★) Sparse な形式的冪級数
- (★★★★) [Bostan-Mori 法](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/bostan_mori.cpp)
- (★★★★) Fiduccia 法 (高速きたまさ法)
- (★★★★) Berlekamp-Massey 法
- (★★★★) 関数の合成
- (★★★★) 逆関数

### 数理最適化 (Optimization)

- (★★☆☆) 二次方程式
- (★★☆☆) 二分探索法 (方程式の解を 1 つ求める)
- (★★☆☆) 三分探索法
- (★★☆☆) 黄金探索法
- (★★★☆) Newton 法
- (★★★★) 単体法 (二段階単体法)
- (★★★★) 分枝限定法


　
<a name="mc"></a>
# 組合せ (MATH : COMBINATORICS)
組合せ論的アルゴリズムたちです

### 二項係数

- (★★☆☆) [二項係数 (オーソドックス, n<=10^7, r<=10^7, p<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient.cpp)
- (★★☆☆) [二項係数 (愚直計算, n<=10^9, r<=10^7, p<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_naive.cpp)
- (★★☆☆) [二項係数 (漸化式計算, n<=5000, r<=5000, p<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_dp.cpp)
- (★★★★) [二項係数 (任意 mod, n<=10^7, r<=10^7, m<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_any_mod.cpp)

### さまざまな数

- (★★☆☆) [重複組合せ](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/combination_with_repetition.cpp)
- (★★★☆) [分割数 P(N, K) (O(NK))](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/partition_number_pnk.cpp)
- (★★★★) [分割数 P(N) (O(N√N))](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/partition_number_pn.cpp)
- (★★★☆) スターリング数
- (★★★☆) ベル数
- (★★★☆) [カタラン数](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/catalan_number.cpp)
- (★★★☆) ベルヌーイ数
- (★★★☆) モンモール数

### 高速なソート

- (★☆☆☆) [クイックソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/quick_sort.cpp)
- (★☆☆☆) [マージソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/merge_sort.cpp)
- (★☆☆☆) [ヒープソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/heap_sort.cpp)
- (★☆☆☆) [計数ソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/counting_sort.cpp)

### さまざまなソート

- (★☆☆☆) [挿入ソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/insertion_sort.cpp)
- (★☆☆☆) [選択ソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/selection_sort.cpp)
- (★☆☆☆) [バブルソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/bubble_sort.cpp)
- (★☆☆☆) [シェルソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/shell_sort.cpp)
- (★☆☆☆) [コムソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/comb_sort.cpp)
- (★☆☆☆) [ボゴソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/bogo_sort.cpp)

### 集合族に関する問題

- (★★★☆) [2-SAT](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/two_sat.cpp)
- (★★★☆) マトロイド上の Greedy 法
- (★★★★) マトロイド交差

### 集合冪級数 (SPS)

- (★★★☆) [高速ゼータ変換](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/fast_zeta_transform.cpp)
- (★★★★) [高速アダマール変換](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/fast_hadmard_transform.cpp)
- (★★★★) AND Convolution
- (★★★★) OR Convolution
- (★★★★) XOR Convolution
- (★★★★) Subset Convolution
- (★★★★) 集合冪級数の exp
- (★★★★) 集合冪級数の合成

### その他

- (★★★☆) [Nim](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/nim.cpp)
- (★★★☆) [LIS and LDS](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/LIS_and_LDS.cpp)
- (★★★☆) [転倒数](https://github.com/drken1215/algorithm/blob/master/DP/inversion_number.cpp)


　
<a name="mmt"></a>
# 整数 (MATH : NUMBER THEORY)
整数論的アルゴリズムたちです

### Modint 

- (★★☆☆) [Modint](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/modint.cpp)
- (★★☆☆) [実行時に法が決まる Modint](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/modint_runtime.cpp)
- (★★★★) [モンゴメリ乗算を用いた Modint (mod は 2^62 未満の奇数)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/modint_montgomery.cpp)

### 約数, 倍数

- (★☆☆☆) [約数列挙](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/divisor.cpp)
- (★☆☆☆) [最大公約数 (Euclid の互除法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/GCD.cpp)
- (★☆☆☆) [最小公倍数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/LCM.cpp)
- (★★☆☆) [拡張 Euclid の互除法](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/extended_GCD.cpp)
- (★☆☆☆) [pn + r (n は非負整数) で表せる整数のうち, x 以上となる最小の整数](https://github.com/drken1215/algorithm/blob/master/Others/amari_lower_bound.cpp)

### 素数

- (★☆☆☆) [素数判定](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/is_prime.cpp)
- (★☆☆☆) [素因数分解](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/prime_factorization.cpp)
- (★★★☆) [Euler のトーティエント関数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/euler_function.cpp)
- (★★★★) [確率的な高速素数判定 (Miller-Rabin 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/is_prime_Miller_Rabin.cpp)
- (★★★★) [確率的な高速素因数分解 (Pollard のロー法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/pollard_rho.cpp)
- (★★★★) [原始根](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/primitive_root.cpp)
- (★★★★) N 以下の素数の個数 (O(N^2/3))

### エラトステネスの篩

- (★☆☆☆) [エラトステネスの篩](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Eratosthenes.cpp)
- (★★☆☆) [エラトステネスの区間篩](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Eratosthenes_segment.cpp)
- (★★★☆) [高速素因数分解, 約数列挙, メビウス関数 (エラトステネスの篩風)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/fast_prime_factorization_eratosthenes.cpp)
- (★★★★) 線形篩
- (★★★★) アトキンの篩

### 方程式

- (★★★☆) [中国剰余定理](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/chinese_reminder_theorem.cpp)
- (★★★☆) [中国剰余定理 (Garner 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/garner.cpp)
- (★★★☆) [離散対数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/mod_log.cpp)
- (★★★★) [ペル方程式](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Pell_equation.cpp)
- (★★★★) 平方剰余 (Tonelli–Shanks 法)

### さまざまな数

- (★★★☆) [有理数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/rational_number.cpp)
- (★★★☆) [多倍長整数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/big_integer.cpp)
- (★★★★) [ガウス整数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/gauss_integer.cpp)

### その他

- (★★☆☆) [平衡三進法展開](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/power_of_three.cpp)
- (★★★☆) [floor sum](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/floor_sum.cpp)
- (★★★☆) [Stern-Brocot 木](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Stern_Brocot.cpp)


　
<a name="st"></a>
# 文字列 (String)
文字列アルゴリズムです

### 構文解析

- (★★★☆) LL(1) 再帰降下パーサ

### 文字列検索

- (★★★☆) [Suffix Array](https://github.com/drken1215/algorithm/blob/master/String/suffix_array.cpp)
- (★★★☆) [ローリングハッシュ](https://github.com/drken1215/algorithm/blob/master/String/rolling_hash.cpp)
- (★★★☆) [単一パターン検索 (KMP 法)](https://github.com/drken1215/algorithm/blob/master/String/knuth_morris_pratt.cpp)
- (★★★☆) 単一パターン検索 (Boyer-Moore 法)
- (★★★★) 複数パターン検索 (Aho-Corasick 法)

### 文字列系アルゴリズム

- (★★★☆) [Z 法](https://github.com/drken1215/algorithm/blob/master/String/z_algorithm.cpp)
- (★★★☆) [Manacher 法](https://github.com/drken1215/algorithm/blob/master/String/manacher.cpp)

### 文字列系データ構造

- (★★★☆) Trie 木
- (★★★★) Palindromic 木 (AOJ 2292)

### その他

- (★★☆☆) [各 index 以降で各文字が最初に登場する index を求める関数](https://github.com/drken1215/algorithm/blob/master/String/next.cpp)


　
<a name="tr"></a>
# 木 (Tree)
木上のクエリに答えるデータ構造や、木に関する問題を解くアルゴリズムの実装です

### 木

- (★★★☆) [木の走査 (部分木サイズ, LCA など)](https://github.com/drken1215/algorithm/blob/master/Tree/run_tree.cpp)
- (★★★☆) [木の直径](https://github.com/drken1215/algorithm/blob/master/Tree/diameter.cpp)
- (★★★☆) 木の重心
- (★★★★) 木の Distance Frequency Table

### 木 DP

- (★★☆☆) 木 DP
- (★★★☆) 全方位木 DP (俗称)
- (★★★☆) 二乗の木 DP (俗称)

### LCA

- (★★★☆) [LCA (ダブリング)](https://github.com/drken1215/algorithm/blob/master/Tree/lca_by_doubling.cpp)
- (★★★☆) [LCA (Euler Tour)](https://github.com/drken1215/algorithm/blob/master/Tree/lca_euler_tour.cpp)
- (★★★☆) [LCA (HL 分解)](https://github.com/drken1215/algorithm/blob/master/Tree/lca_heavy_light_decomposition.cpp)

### 木上のクエリ処理

- (★★★☆) [Euler Tour (頂点上のクエリ)](https://github.com/drken1215/algorithm/blob/master/Tree/euler_tour_on_nodes.cpp)
- (★★★☆) [Euler Tour (辺上のクエリ)](https://github.com/drken1215/algorithm/blob/master/Tree/euler_tour_on_edges.cpp)
- (★★★☆) [HL 分解](https://github.com/drken1215/algorithm/blob/master/Tree/heavy_light_decomposition.cpp)
- (★★★☆) [重心分解](https://github.com/drken1215/algorithm/blob/master/Tree/tree_centroid_decomposition.cpp)
- (★★★☆) DSU on Tree
- (★★★★) Link-Cut 木
- (★★★★) toptree

### その他の問題

- (★★★☆) [強平衡二分木の Distance Frequency Table](https://github.com/drken1215/algorithm/blob/master/Tree/find_various_values_of_binary_tree.cpp)
- Level Ancester


　
<a name="ot"></a>
# その他 (OTHERS)
その他のアルゴリズムです

### 入出力

- (★☆☆☆) [デバッグストリーム, chmin, chmax](https://github.com/drken1215/algorithm/blob/master/Others/debug.cpp)
- (★★★☆) Nyaan's 高速入出力

### グリッド

- (★☆☆☆) [グリッドの 4 近傍, 8 近傍](https://github.com/drken1215/algorithm/blob/master/Others/grid_neighbors.cpp)
- (★☆☆☆) [ハニカムの 6 近傍](https://github.com/drken1215/algorithm/blob/master/Others/honeycomb_neighbors.cpp)

### ビット演算

- (★★☆☆) [XorShift, ランダムシャッフル](https://github.com/drken1215/algorithm/blob/master/Others/xorshift.cpp)
- (★★☆☆) [next_combination](https://github.com/drken1215/algorithm/blob/master/Others/next_combination.cpp)
- (★★☆☆) [部分集合の部分集合](https://github.com/drken1215/algorithm/blob/master/Others/subset_enumeration.cpp)

### 探索法

- (★★★☆) α-β 探索
- (★★★☆) 焼き鈍し法
- (★★★☆) A*
- (★★★☆) IDA*
- (★★★☆) Baby-Step Giant-Step 法

### その他

- (★★☆☆) [タイマー](https://github.com/drken1215/algorithm/blob/master/Others/timer.cpp)
- (★★☆☆) [サイコロ](https://github.com/drken1215/algorithm/blob/master/Others/dice.cpp)
- (★★☆☆) [曜日](https://github.com/drken1215/algorithm/blob/master/Others/day_of_the_week.cpp)
- (★★★☆) 四面体 (AOJ 2060)


　
# License
These codes are licensed under CC0.
[![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png "CC0")](http://creativecommons.org/publicdomain/zero/1.0/deed.ja)

