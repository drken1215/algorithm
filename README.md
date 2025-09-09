# 様々なアルゴリズムの実装例
データ構造や数論的アルゴリズムまで、様々な分野のアルゴリズムたちを C++23 で実装しています。  
アルゴリズム系の研究開発において計算機実験が必要になる場面や、
プログラミングコンテストに参加する場面などを想定して、
「実装例」または「ライブラリ」として使用することを念頭に置いています。

　

|分類|内容|具体例|
|---|---|---|
|**[DATA STRUCTURE](#ds)**|各種データ構造|Union-Find、Sparse Table など|
|**[DATA STRUCTURE : SEGMENT](#dss)**|区間クエリに強いデータ構造|セグメント木、BIT など|
|**[GEOMETRY](#ge)**|計算幾何|円の交点など|
|**[GRAPH](#gt)**|グラフアルゴリズム|強連結成分分解など|
|**[GRAPH : NETWORK FLOW](#gnf)**|ネットワークフローアルゴリズム|Ford-Fulkerson 法など|
|**[MATH : ALGEBRA](#ma)**|代数的アルゴリズム|行列計算など|
|**[MATH : COMBINATORICS](#mc)**|組合せ論的アルゴリズム|modint、Nim など|
|**[MATH : NUMBER THEORY](#mmt)**|整数論的アルゴリズム|素因数分解、最大公約数など|
|**[OPTIMIZATION](#opt)**|最適化や探索のアルゴリズム|二分探索, 三分探索など|
|**[STRING](#st)**|文字列アルゴリズム|Suffix Array、KMP 法など|
|**[TREE](#tr)**|木上のデータ構造とアルゴリズム|Euler ツアー、木の直径など|
|**[OTHERS](#ot)**|その他|xorshift、サイコロなど|

## 難易度表記の目安

- (★☆☆☆)：一般教養、NoviSteps グレード基準で 2Q 以下
- (★★☆☆)：初等典型、NoviSteps グレード基準で 1Q, 1D
- (★★★☆)：中堅典型、NoviSteps グレード基準で 2D, 3D
- (★★★★)：高度典型、NoviSteps グレード基準で 4D 以上


　

<a name="ds"></a>

# データ構造 (DATA STRUCTURE)
各種データ構造の実装です

## Union-Find

- (★☆☆☆) [Union-Find (union by size)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_size.cpp)
- (★☆☆☆) [Union-Find (union by rank)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_rank.cpp)
- (★★☆☆) [重みつき Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/weighted_union_find.cpp)
- (★★★☆) [列挙可能 Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_enumerable.cpp)
- (★★★☆) [undo つき Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_can_undo.cpp)
- (★★★☆) [部分永続 Union-Find](https://github.com/drken1215/algorithm/blob/master/DataStructure/partially_persistent_union_find.cpp)

## キュー, ヒープ

- (★☆☆☆) [二分ヒープ](https://github.com/drken1215/algorithm/blob/master/DataStructure/heap.cpp)
- (★★☆☆) [削除可能 priority queue](https://github.com/drken1215/algorithm/blob/master/DataStructure/removable_queue.cpp)
- (★★★☆) [両端 priority queue (削除も可能)](https://github.com/drken1215/algorithm/blob/master/DataStructure/double_ended_priority_queue.cpp)
- (★★★★) Skew Heap (マージ可能ヒープ)
- (★★★★) Paring Heap (マージ可能ヒープ)
- (★★★★) Radix Heap
- (★★★★) Fibonacci Heap

## ハッシュ

- (★★☆☆) [ハッシュ](https://github.com/drken1215/algorithm/blob/master/DataStructure/hash.cpp)
- (★★☆☆) [Zobrist ハッシュ](https://github.com/drken1215/algorithm/blob/master/DataStructure/zobrist_hash.cpp)
- (★★★☆) 根付き木のハッシュ

## N 以下の非負整数値の順序つき集合

- (★★☆☆) [BIT と BIT 上の二分探索](https://github.com/drken1215/algorithm/blob/master/DataStructure/predecessor_bit.cpp)
- (★★★★) [64 分木](https://github.com/drken1215/algorithm/blob/master/DataStructure/predecessor_64_tree.cpp)
- (★★★★) van Emde Boas 木

## 永続データ構造

- (★★★★) 永続配列
- (★★★★) 完全永続 Union-Find
- (★★★★) 永続キュー
- (★★★★) 永続セグメント木
- (★★★★) 永続赤黒木

## その他

- (★★★☆) [可変長 BitSet](https://github.com/drken1215/algorithm/blob/master/DataStructure/advanced_bitset.cpp)
- (★★★☆) [並列二分探索](https://github.com/drken1215/algorithm/blob/master/DataStructure/parallel_binary_search.cpp)
- (★★★☆) [Cartesian 木](https://github.com/drken1215/algorithm/blob/master/DataStructure/cartesian_tree.cpp)
- (★★★★) [Binary Trie](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_trie.cpp)
- (★★★★) Dynamic Connectivity


　

<a name="dss"></a>

# 区間系データ構造 (DATA STRUCTURE : SEGMENT)
セグメント木や BIT など、区間クエリに強いデータ構造の実装です

## セグメント木, Binary Indexed 木

- (★★☆☆) [セグメント木](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/segment_tree.cpp)
- (★★★☆) [セグメント木 (遅延評価)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/segment_tree_lazy.cpp)
- (★★☆☆) [BIT](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_indexed_tree.cpp)
- (★★★☆) [BIT (区間加算, 区間和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_indexed_tree_RAQ.cpp)
- (★★★★) Segment Tree Beats (俗称)

## セグメント木の応用
- (★★☆☆) [RMQ (セグメント木)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/range_minimum_query.cpp)
- (★★★☆) [Starry Sky 木 (俗称)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/starry_sky_tree.cpp)
- (★★★☆) [マージソート過程木 (領域木)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/merge_sort_tree.cpp)

## 動的・二次元セグメント木

- (★★★☆) 動的セグメント木
- (★★★☆) [二次元セグメント木](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/segment_tree_2d.cpp)
- (★★★☆) [二次元 BIT](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_indexed_tree_2D.cpp)
- (★★★★) [二次元 BIT (領域加算, 領域和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/binary_indexed_tree_2D_RAQ.cpp)
- (★★★★) 動的二次元セグメント木
- (★★★★) 動的二次元 BIT

## Sparse Table

- (★★★☆) [Sparse Table](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/sparse_table.cpp)
- (★★★★) [Disjoint Sparse Table](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/disjoint_sparse_table.cpp)
- (★★★★) 二次元 Sparse Table

## ウェーブレット行列

- (★★★★) [ウェーブレット行列](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/wavelet_matrix.cpp)
- (★★★★) [BIT on ウェーブレット行列 (一点加算対応)](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/BIT_on_wavelet_matrix.cpp)
- (★★★★) [セグメント木 on ウェーブレット行列](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/segment_tree_on_wavelet_matrix.cpp)
- (★★★★) 動的ウェーブレット行列

## 平衡二分探索木

- (★★★★) [RBST](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/randomized_binary_search_tree.cpp)
- (★★★★) Treap 木
- (★★★★) AVL 木
- (★★★★) Splay 木
- (★★★★) 赤黒木

## 各種高速化アルゴリズム

- (★☆☆☆) [累積和](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/cumulative_sum.cpp)
- (★☆☆☆) [二次元累積和](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/cumulative_sum_2D.cpp)
- (★★★☆) [スライド最小値](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/sliding_minimum.cpp)
- (★★★☆) [SWAG](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/swag.cpp)
- (★★★☆) [Mo 法](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/mo.cpp)

## その他

- (★★★☆) [区間の集合を set で管理する（del and add must be invertible）](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/intervals_management.cpp)
- (★★★☆) [区間の集合を set で管理する（del and add can be non-invertible）](https://github.com/drken1215/algorithm/blob/master/DataStructureSegment/intervals_management_with_noninvertible_del.cpp)





　
<a name="ge"></a>
# 幾何 (GEOMETRY)
幾何ライブラリです

- (★★★☆) [全部乗せ](https://github.com/drken1215/algorithm/blob/master/Geometry/All.cpp)

## 基本要素

- (★☆☆☆) [基本要素 (点, 線分, 円)](https://github.com/drken1215/algorithm/blob/master/Geometry/basic_elements.cpp)
- (★☆☆☆) [偏角ソート](https://github.com/drken1215/algorithm/blob/master/Geometry/arg_sort.cpp)

## 点, 線分, 三角形などの位置関係

- (★☆☆☆) [点と線分の位置関係 (ccw)](https://github.com/drken1215/algorithm/blob/master/Geometry/ccw.cpp)
- (★☆☆☆) [点と三角形の包含関係](https://github.com/drken1215/algorithm/blob/master/Geometry/is_contain_in_the_triangle.cpp)

## 射影, 交差判定, 距離

- (★☆☆☆) [射影](https://github.com/drken1215/algorithm/blob/master/Geometry/projection.cpp)
- (★☆☆☆) [線分と線分の交差判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_intersect_two_segments.cpp)
- (★☆☆☆) [線分と線分との距離](https://github.com/drken1215/algorithm/blob/master/Geometry/distance_two_segments.cpp)

## 直線や円の交点

- (★☆☆☆) [直線と直線の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_two_lines.cpp)
- (★★☆☆) [円と直線の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_line_circle.cpp)
- (★★☆☆) [円と円の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_two_circles.cpp)
- (★★☆☆) [円と線分の交点](https://github.com/drken1215/algorithm/blob/master/Geometry/crosspoint_segment_circle.cpp)

## 接線

- (★★☆☆) [接線](https://github.com/drken1215/algorithm/blob/master/Geometry/tanline.cpp)
- (★★★☆) [共通接線 (2 円)](https://github.com/drken1215/algorithm/blob/master/Geometry/common_tanline.cpp)

## 多角形

- (★☆☆☆) [多角形の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_polygon.cpp)
- (★★☆☆) [点と多角形の包含判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_contain_in_the_polygon.cpp)
- (★★☆☆) [凸性判定](https://github.com/drken1215/algorithm/blob/master/Geometry/is_convex.cpp)
- (★★☆☆) [凸包](https://github.com/drken1215/algorithm/blob/master/Geometry/convex_hull.cpp)
- (★★★☆) [凸多角形の直径](https://github.com/drken1215/algorithm/blob/master/Geometry/diameter.cpp)
- (★★★☆) [凸多角形の切断](https://github.com/drken1215/algorithm/blob/master/Geometry/convex_cut.cpp)
- (★★★☆) 凸多角形と直線の交点 O(log N)
- (★★★☆) [円と円の共通部分の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_common_two_circles.cpp)
- (★★★☆) [円と多角形との共通部分の面積](https://github.com/drken1215/algorithm/blob/master/Geometry/area_common_circle_polygon.cpp)
- (★★★☆) [ボロノイ図 (単純ver, O(N^3))](https://github.com/drken1215/algorithm/blob/master/Geometry/voronoi.cpp)
- (★★★★) ボロノイ図 (単純ver, O(N log N))
- (★★★★) ドロネーの三角形分割 (期待値 O(N^2))

## 三次元幾何

- (★★★☆) [三次元幾何一式](https://github.com/drken1215/algorithm/blob/master/Geometry/basic_elements_3D.cpp)

## その他

- (★★☆☆) [垂直二等分線](https://github.com/drken1215/algorithm/blob/master/Geometry/perpendicular_bisector.cpp)
- (★★★☆) [3 点を通る円](https://github.com/drken1215/algorithm/blob/master/Geometry/circumscribed_circle.cpp)
- (★★★☆) [アポロニウスの円](https://github.com/drken1215/algorithm/blob/master/Geometry/apollonius.cpp)
- (★★★☆) [最近点対](https://github.com/drken1215/algorithm/blob/master/Geometry/closest_two_points.cpp)
- (★★★★) 最小包含円
- (★★★☆) 線分併合
- (★★★☆) 線分アレンジメント
- (★★★☆) 双対変換
- (★★★★) kd 木


　

<a name="gt"></a>
# グラフ (GRAPH)
グラフアルゴリズムです

## DFS, BFS

- (★☆☆☆) [連結成分の個数 (DFS)](https://github.com/drken1215/algorithm/blob/master/Graph/dfs.cpp)
- (★☆☆☆) [連結成分の個数 (BFS)](https://github.com/drken1215/algorithm/blob/master/Graph/bfs.cpp)
- (★☆☆☆) [二部グラフ判定 (DFS)](https://github.com/drken1215/algorithm/blob/master/Graph/is_bipartite_dfs.cpp)
- (★☆☆☆) [二部グラフ判定 (BFS)](https://github.com/drken1215/algorithm/blob/master/Graph/is_bipartite_bfs.cpp)
- (★★☆☆) [トポロジカルソート (DFS)](https://github.com/drken1215/algorithm/blob/master/Graph/topological_sort_dfs.cpp)
- (★★☆☆) [トポロジカルソート (BFS)](https://github.com/drken1215/algorithm/blob/master/Graph/topological_sort_bfs.cpp)

## 連結成分分解

- (★☆☆☆) [閉路検出 (サイクル検出)](https://github.com/drken1215/algorithm/blob/master/Graph/cycle_detection.cpp)
- (★☆☆☆) [Functional グラフの閉路検出 (サイクル検出)](https://github.com/drken1215/algorithm/blob/master/Graph/functional_graph_cycle_detection.cpp)
- (★★☆☆) [強連結成分分解](https://github.com/drken1215/algorithm/blob/master/Graph/strongly_connected_components.cpp)
- (★★★☆) [橋, 関節点列挙 (Low-Link)](https://github.com/drken1215/algorithm/blob/master/Graph/low_link.cpp)
- (★★★☆) [二重辺連結成分分解 (Bridge-Block 木)](https://github.com/drken1215/algorithm/blob/master/Graph/two_edge_connected_components.cpp)
- (★★★☆) [二重頂点連結成分分解 (Block-Cut 木)](https://github.com/drken1215/algorithm/blob/master/Graph/biconnected_components.cpp)
- (★★★★) 三重辺連結成分分解 (SPQR 木)

## 最短路問題 (基本)

- (★☆☆☆) [重みなしグラフの最短路 (BFS, O(E))](https://github.com/drken1215/algorithm/blob/master/Graph/shortest_path_bfs.cpp)
- (★☆☆☆) [重みが 0, 1 のみのグラフの最短路 (0-1 BFS, O(E))](https://github.com/drken1215/algorithm/blob/master/Graph/shortest_path_01bfs.cpp)
- (★☆☆☆) [単一始点最短路 (Dijkstra 法, 正辺のみ, O(V + E log V))](https://github.com/drken1215/algorithm/blob/master/Graph/shortest_path_dijkstra.cpp)
- (★☆☆☆) [単一始点最短路 (Bellman-Ford 法, 負辺対応, O(VE))](https://github.com/drken1215/algorithm/blob/master/Graph/shortest_path_bellman_ford.cpp)
- (★☆☆☆) [全頂点対間最短路 (Floyd-Warshall 法, O(V^3))](https://github.com/drken1215/algorithm/blob/master/Graph/floyd_warshall.cpp)
- (★★★☆) 全頂点対間最短路 (Johnson 法, O(EV log V))
- (★★★☆) SPFA
- (★★★★) 補グラフの最短路

## 全域木, 路に関する問題

- (★★☆☆) 最小全域木 (Kruskal 法)
- (★★★☆) 有向 Euler 路
- (★★★☆) [無向 Euler 路](https://github.com/drken1215/algorithm/blob/master/Graph/euler_tour_undirected.cpp)
- (★★★★) 最小有向全域木 (Chu-Liu/Edmonds 法)
- (★★★★) [無向グラフの全域木の個数 (行列木定理)](https://github.com/drken1215/algorithm/blob/master/Graph/matrix_tree_theorem.cpp)
- (★★★★) 最小シュタイナー木 (O(V 3^t + V^2 2^t + V^3))

## グラフ上の有名問題

- (★★★☆) [最大安定集合問題 (O(1.381^V))](https://github.com/drken1215/algorithm/blob/master/Graph/maximum_stable_set.cpp)
- (★★★☆) 最大クリーク列挙（O(1.443^V)）
- (★★★☆) [頂点彩色 (O(2^V V))](https://github.com/drken1215/algorithm/blob/master/Graph/vertex_coloring.cpp)
- (★★★★) 辺彩色
- (★★★★) 二部グラフの辺彩色 (Alon 法, O(E log E))


　
　
<a name="gnf"></a>
# ネットワークフロー (GRAPH : NETWORK FLOW)
グラフネットワークフロー関連のアルゴリズムです

## 最大流

- (★★★☆) [最大流 (Ford-Fulkerson 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/max_flow_ford_fulkerson.cpp)
- (★★★☆) [最大流 (Dinic 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/max_flow_dinic.cpp)

## 最小費用流

- (★★★☆) [最小費用流 (Primal-Dual 法, 正辺のみ)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_flow_primal_dual.cpp)
- (★★★☆) [最小費用流 (Primal-Dual 法, 負辺対応)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_flow_primal_dual_negative.cpp)
- (★★★★) [最小費用循環流 (Goldberg-Tarjan 法, by cost-scaling, 負閉路 OK)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_circulating_flow.cpp)
- (★★★★) [最小費用 b-flow (by cost-scaling)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_b_flow_by_cost_scaling.cpp)
- (★★★★) [最小費用 b-flow (by ネットワーク単体法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cost_b_flow_by_network_simplex_method.cpp)
- (★★★★) 最小凸費用流

## 劣モジュラ関数のグラフ表現

- (★★★☆) [Project Selection Problem (俗称: 燃やす埋める)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/project_selection_problem.cpp)
- (★★★☆) [2 変数劣モジュラ関数の和を表す最小カット (燃やす埋めるの一般化)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/two_variable_submodular_optimization.cpp)
- (★★★★) [3 変数劣モジュラ関数の和を表す最小カット (燃やす埋めるの一般化)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/three_variable_submodular_optimization.cpp)

## カット

- (★★★☆) [最小カット (= 最大流)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/min_cut.cpp)
- (★★★★) 全域最小カット（Stoer-Wanger 法）
- (★★★★) 全頂点対間最小カット (Nagamochi-Ibaraki 法)
- (★★★★) Gomory-Hu 木

## マッチング

- (★★★☆) [二部マッチング (Hopcroft-Karp 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/hopcroft_karp.cpp)
- (★★★☆) 重みつき二部マッチング (Hungarian 法)
- (★★★★) 一般グラフの最大マッチング (Edmonds 法)
- (★★★★) 一般グラフの最大マッチング (行列補間)
- (★★★★) 重み付き一般グラフの最大マッチング


　

<a name="ma"></a>
# 代数 (MATH : ALGEBRA)
行列計算など代数的計算に関するアルゴリズムです

## 行列

- (★★★☆) [行列累乗, 連立一次方程式 (実数)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_double.cpp)
- (★★★☆) [行列累乗, 行列式, 逆行列, ランク, 連立一次方程式 (mod. p)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_modp.cpp)
- (★★★☆) [行列累乗, 行列式, 逆行列, ランク, 連立一次方程式 (binary)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_binary.cpp)
- (★★★☆) [任意 mod 行列式 (by Euclid の互助法を用いた掃き出し法)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix_determinant_in_general_mod.cpp)
- (★★★★) [特性多項式の計算 (O(N^3))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/characteristic_polynomial.cpp)
- (★★★★) ハフニアン (完全マッチングの個数に帰着)

## 多項式行列

- (★★★★) 多項式行列の行列式 (by 除算を全く行わない方法で, O(N^4))
- (★★★★) 多項式行列の行列式 (by ラグランジュ補間によって, O(N^3 d))
- (★★★★) 多項式行列の行列式 (by Euclid の互助法を用いた掃き出し法)
- (★★★★) [行列式 det(M0 + M1x) の計算 (O(N^3))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/determinant_matrix_linear_expression.cpp)
- (★★★★) 多項式行列の prefix product M(0)M(1)...M(K-1)

## 特殊な行列の乗法、逆行列

- (★★★★) Black Box Linear Algebra (行列式などが O(N^2 + N T(N)))
- (★★★★) 巡回行列 (行列式などが O(N^2))
- (★★★★) 上三角 Toeplitz 行列 (行列式などが O(N^2))
- (★★★★) K 重対角行列 (行列式などが O(NK^2))

## FFT, NTT, Convolution

- (★★★☆) [FFT (高速フーリエ変換)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/FFT.cpp)
- (★★★☆) [NTT (高速剰余変換)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/NTT.cpp)
- (★★★☆) [任意 mod 畳み込み (mod < 10^9)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/NTT_any_mod.cpp)
- (★★★★) [mod 2^64 畳み込み](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/NTT_mod_2_64.cpp)

## 形式的冪級数 (FPS)

- (★★★★) [形式的冪級数：全部乗せ](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/formal_power_series.cpp)
- (★★★★) [Inv of FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/inv_formal_power_series.cpp)
- (★★★★) [Exp of FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/exp_formal_power_series.cpp)
- (★★★★) [Log of FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/log_formal_power_series.cpp)
- (★★★★) [Pow of FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/pow_formal_power_series.cpp)
- (★★★★) [Sqrt of FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/sqrt_formal_power_series.cpp)

## Sparse な形式的冪級数 (FPS)

- (★★★★) [Inv of Sparse FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/inv_sparse_formal_power_series.cpp)
- (★★★★) [Exp of Sparse FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/exp_sparse_formal_power_series.cpp)
- (★★★★) [Log of Sparse FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/log_sparse_formal_power_series.cpp)
- (★★★★) [Pow of Sparse FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/pow_sparse_formal_power_series.cpp)
- (★★★★) [Sart of Sparse FPS](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/sqrt_sparse_formal_power_series.cpp)

## FPS のアルゴリズム

- (★★★★) [Bostan-Mori 法](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/bostan_mori.cpp)
- (★★★★) [Power Projection](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/power_projection.cpp)
- (★★★★) [FPS の合成 (Kinoshita-Li 法, O(N (log N)^2)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/composition_formal_power_series.cpp)
- (★★★★) [FPS の逆関数 (Kinoshita-Li 法, O(N (log N)^2)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/compositional_inverse_formal_power_series.cpp)
- (★★★★) Berlekamp-Massey 法

## 多項式の基底変換

- (★★★☆) [Polynomial Taylor Shift (O(N log N))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_taylor_shift.cpp)
- (★★★★) [Multipoint Evaluation (O(M(log M)^2 + N log N))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/multipoint_evaluation.cpp)
- (★★★★) [Multipoint Evaluation (等比数列のとき) (by chirp z-transform, O(N log N))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/multipoint_evaluation_in_geometric_sequence.cpp)
- (★★★★) [多項式補間 (O(N(log N)^2))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_interpolation.cpp)
- (★★★★) [多項式補間 (等比数列のとき) (O(N log N))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_interpolation_in_geometric_sequence.cpp)
- (★★★★) [多項式の基底変換：Monomial 基底 → Newton 基底 (O(N(log N)^2))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_monomial_to_newton.cpp)

## 多項式のその他のアルゴリズム

- (★★★☆) [多項式の累乗 f(x)^e mod g(x)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_mod_pow.cpp)
- (★★★☆) [多項式の middle product (c[i] = sum_j a[i+j]b[j])](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_middle_product.cpp)
- (★★★★) [多項式の評価点シフト (O((N + M)log(N + M))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_shift_sampling.cpp)
- (★★★★) [多項式の除算 (by NTT, inv of FPS, O(N log N))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_div.cpp)
- (★★★★) [多項式 GCD (by half-gcd) (O(N(log N)^2))](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_gcd.cpp)
- (★★★★) [多項式の零点を求める (mod 998244353)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/polynomial_root_finding.cpp)

## さまざまな和

- (★★★☆) [floor sum](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/floor_sum.cpp)
- (★★★★) Bk = Σ{i=0}^{n-1} i^k を係数とする多項式
- (★★★★) Σ{i=0}^{n-1} r^i i^d
- (★★★★) Σ{i=0}^{∞} r^i i^d
- (★★★★) Karatsuba 法



　
<a name="mc"></a>
# 組合せ (MATH : COMBINATORICS)
組合せ論的アルゴリズムたちです

## 二項係数

- (★☆☆☆) [二項係数 (オーソドックス, n<=10^7, r<=10^7, p<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient.cpp)
- (★☆☆☆) [二項係数 (愚直計算, n<=10^9, r<=10^7, p<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_naive.cpp)
- (★☆☆☆) [二項係数 (漸化式計算, n<=5000, r<=5000, p<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_dp.cpp)
- (★★★☆) [二項係数 (任意 mod, n<=10^7, r<=10^7, m<=10^9)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_binomial_coefficient_any_mod.cpp)

## さまざまな数

- (★★☆☆) [重複組合せ](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/combination_with_repetition.cpp)
- (★★★☆) [カタラン数](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/catalan_number.cpp)
- (★★★☆) [分割数 P(N, K) (O(NK))](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/partition_number_pnk.cpp)
- (★★★★) [分割数 P(N) (O(N√N))](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/partition_number_pn.cpp)
- (★★★☆) スターリング数
- (★★★☆) ベル数
- (★★★☆) ベルヌーイ数
- (★★★☆) モンモール数

## 高速なソート

- (★☆☆☆) [クイックソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/quick_sort.cpp)
- (★☆☆☆) [マージソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/merge_sort.cpp)
- (★☆☆☆) [ヒープソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/heap_sort.cpp)
- (★☆☆☆) [計数ソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/counting_sort.cpp)

## さまざまなソート

- (★☆☆☆) [挿入ソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/insertion_sort.cpp)
- (★☆☆☆) [選択ソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/selection_sort.cpp)
- (★☆☆☆) [バブルソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/bubble_sort.cpp)
- (★☆☆☆) [シェルソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/shell_sort.cpp)
- (★☆☆☆) [コムソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/comb_sort.cpp)
- (★☆☆☆) [ボゴソート](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/bogo_sort.cpp)

## 集合族に関する問題

- (★★★☆) [2-SAT](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/two_sat.cpp)
- (★★★☆) マトロイド上の Greedy 法
- (★★★★) マトロイド交差

## 集合冪級数 (SPS)

- (★★★☆) [高速ゼータ変換](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/fast_zeta_transform.cpp)
- (★★★★) [高速アダマール変換](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/fast_hadmard_transform.cpp)
- (★★★★) AND Convolution
- (★★★★) OR Convolution
- (★★★★) XOR Convolution
- (★★★★) Subset Convolution
- (★★★★) 集合冪級数の exp
- (★★★★) 集合冪級数の合成

## その他

- (★★☆☆) [Nim](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/nim.cpp)
- (★★☆☆) [LIS and LDS](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/LIS_and_LDS.cpp)
- (★★☆☆) [転倒数](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/inversion_number.cpp)


　
<a name="mmt"></a>
# 整数 (MATH : NUMBER THEORY)
整数論的アルゴリズムたちです

## Mod の処理

- (★☆☆☆) [a^n, a^{-1} (mod m)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/mod_pow_inv.cpp)
- (★☆☆☆) [Modint](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/modint.cpp)
- (★★☆☆) [実行時に法が決まる Modint](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/modint_runtime.cpp)
- (★★★☆) [モンゴメリ乗算を用いた Modint (mod は 2^62 未満の奇数)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/modint_montgomery.cpp)

## 約数, 倍数

- (★☆☆☆) [約数列挙](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/divisor.cpp)
- (★☆☆☆) [最大公約数 (Euclid の互除法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/GCD.cpp)
- (★☆☆☆) [最小公倍数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/LCM.cpp)
- (★☆☆☆) [拡張 Euclid の互除法](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/extended_GCD.cpp)

## 素数

- (★☆☆☆) [素数判定](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/is_prime.cpp)
- (★☆☆☆) [素因数分解](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/prime_factorization.cpp)
- (★★☆☆) [Euler のトーティエント関数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/euler_function.cpp)
- (★★★☆) [確率的な高速素数判定 (Miller-Rabin 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/is_prime_Miller_Rabin.cpp)
- (★★★☆) [確率的な高速素因数分解 (Pollard のロー法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/pollard_rho.cpp)
- (★★★☆) [原始根](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/primitive_root.cpp)
- (★★★★) [位数 (a^x ≡ 1 (mod p) となる最小の x)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/order.cpp)
- (★★★★) N 以下の素数の個数 (O(N^2/3))

## エラトステネスの篩、GCD Convolution

- (★☆☆☆) [エラトステネスの篩](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Eratosthenes.cpp)
- (★★☆☆) [エラトステネスの区間篩](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Eratosthenes_segment.cpp)
- (★★☆☆) [高速素因数分解, 約数列挙, メビウス関数 (エラトステネスの篩風)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/fast_prime_factorization_eratosthenes.cpp)
- (★★★★) [添字 GCD 畳み込み](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/fast_gcd_convolution.cpp)
- (★★★★) 線形篩
- (★★★★) アトキンの篩

## 方程式

- (★★★☆) [中国剰余定理](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/chinese_reminder_theorem.cpp)
- (★★★☆) [中国剰余定理 (Garner 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/garner.cpp)
- (★★★☆) [離散対数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/mod_log.cpp)
- (★★★★) [ペル方程式](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Pell_equation.cpp)
- (★★★★) [平方剰余 (Tonelli–Shanks 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/mod_sqrt.cpp)
- (★★★★) [Kth Root (MOD) (Tonelli–Shanks 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/mod_kth_root.cpp)

## さまざまな数

- (★★☆☆) [128 ビット整数 (ラッパー)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/int128.cpp)
- (★★★☆) [多倍長整数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/big_integer.cpp)
- (★★★☆) [有理数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/rational_number.cpp)
- (★★★★) [ガウス整数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/gauss_integer.cpp)

## その他

- (★☆☆☆) [m で割って r 余る, x 以上の最小の整数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/amari_lower_bound.cpp)
- (★☆☆☆) [平衡三進法展開](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/power_of_three.cpp)
- (★★☆☆) [x^K <= N となる最大の整数 x](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/kth_root.cpp)
- (★★★☆) [Stern-Brocot 木](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/Stern_Brocot.cpp)


　

<a name="opt"></a>
# 最適化, 探索 (OPTIMIZATION, SEARCH)
最適化や探索に関するアルゴリズムです

## さまざまな全探索

- (★★☆☆) [next_combination (nCk 通りの全探索)](https://github.com/drken1215/algorithm/blob/master/Optimization/next_combination.cpp)
- (★★☆☆) [部分集合の部分集合 (3^N 通りの全探索)](https://github.com/drken1215/algorithm/blob/master/Optimization/subset_enumeration.cpp)
- (★★☆☆) [数独ソルバー](https://github.com/drken1215/algorithm/blob/master/Optimization/sudoku.cpp)

## 動的計画法

- (★☆☆☆) [ナップサック問題](https://github.com/drken1215/algorithm/blob/master/Optimization/knapsack.cpp)
- (★☆☆☆) [LIS](https://github.com/drken1215/algorithm/blob/master/Optimization/longest_increasing_sequence.cpp)
- (★☆☆☆) [LCS](https://github.com/drken1215/algorithm/blob/master/Optimization/lcs.cpp)
- (★☆☆☆) [編集距離](https://github.com/drken1215/algorithm/blob/master/Optimization/edit_distance.cpp)
- (★★☆☆) [グリッドに含まれる最大正方形](https://github.com/drken1215/algorithm/blob/master/Optimization/largest_square_in_grid.cpp)
- (★★★☆) [ヒストグラム長方形面積最大化](https://github.com/drken1215/algorithm/blob/master/Optimization/histogram.cpp)
- (★★★★) 最適二分探索木 (O(N^2), Monge 性を活かした区間 DP)
- (★★★★) 最適二分探索木 (O(N log N), Hu-Tucker 法)

## Convex Hull Trick

- (★★★☆) [Convex Hull Trick (直線：傾き単調, O(log N)) (クエリも単調なら O(1))](https://github.com/drken1215/algorithm/blob/master/Optimization/convex_hull_trick_slope_monotone.cpp)
- (★★★★) [Convex Hull Trick (直線：傾き単調でなくてよい, Li Chao Tree, O(log N))](https://github.com/drken1215/algorithm/blob/master/Optimization/convex_hull_trick.cpp)
- (★★★★) [一般化 Convex Hull Trick (Monge, O(log N))](https://github.com/drken1215/algorithm/blob/master/Optimization/convex_hull_trick_general.cpp)

## Monge 性を活用する DP 高速化技法

- (★★★☆) [Monotone 行最小値問題 (Monotone Minima, O(H + W log H))](https://github.com/drken1215/algorithm/blob/master/Optimization/monotone_minima.cpp)
- (★★★☆) [Monotone 単一始点最短路問題 (D&D Monotone Minima, O(N (log N)^2))](https://github.com/drken1215/algorithm/blob/master/Optimization/monotone_minima.cpp)
- (★★★☆) TM 行最小値問題 (SMAWK, O(H + W))
- (★★★☆) TM 単一始点最短路問題 (D&D SMAWK, O(N log N))
- (★★★★) Monge 単一始点最短路問題（by noshi's 簡易 LARSCH, O(N log N))
- (★★★★) Monge 単一始点最短路問題（by LARSCH, O(N))
- (★★★★) Monge グラフ上の d-辺最短路
- (★★★★) Aliens DP

## その他の DP 高速化技法

- (★★★★) Slope Trick

## 数理最適化

- (★☆☆☆) 二次方程式
- (★☆☆☆) 二分探索法 (方程式の解を 1 つ求める)
- (★★☆☆) 三分探索法
- (★★☆☆) 黄金探索法
- (★★★☆) Newton 法
- (★★★★) [単体法 (二段階単体法)](https://github.com/drken1215/algorithm/blob/master/Optimization/simplex_method.cpp)
- (★★★★) 分枝限定法

## さまざまな探索法

- (★★★☆) α-β 探索
- (★★★☆) 焼き鈍し法
- (★★★☆) A*
- (★★★☆) IDA*

## 指数時間アルゴリズム

- (★★★★) Set Cover
- (★★★★) k-Cover (O(2^N N))
- (★★★★) k-partition (O(2^N N^3))


　

<a name="st"></a>
# 文字列 (String)
文字列アルゴリズムです

### 構文解析

- (★★☆☆) [LL(1) 再帰下降パーサ](https://github.com/drken1215/algorithm/blob/master/String/parser.cpp)

### 文字列検索

- (★★☆☆) [ローリングハッシュ](https://github.com/drken1215/algorithm/blob/master/String/rolling_hash.cpp)
- (★★★☆) [Suffix Array](https://github.com/drken1215/algorithm/blob/master/String/suffix_array.cpp)
- (★★★☆) [単一パターン検索 (KMP 法)](https://github.com/drken1215/algorithm/blob/master/String/knuth_morris_pratt.cpp)
- (★★★☆) 単一パターン検索 (Boyer-Moore 法)
- (★★★★) 複数パターン検索 (Aho-Corasick 法)

### さまざまな文字列アルゴリズム

- (★★☆☆) [Z 法](https://github.com/drken1215/algorithm/blob/master/String/z_algorithm.cpp)
- (★★★☆) [Manacher 法](https://github.com/drken1215/algorithm/blob/master/String/manacher.cpp)

### さまざまな文字列データ構造

- (★★★☆) Trie 木
- (★★★★) Palindromic 木 (AOJ 2292)

### その他

- (★★☆☆) [各 index 以降で各文字が最初に登場する index を求める関数](https://github.com/drken1215/algorithm/blob/master/String/next.cpp)
- (★★★☆) [文字列の相異なる部分列の個数 (O(cN))](https://github.com/drken1215/algorithm/blob/master/String/num_of_subsequences.cpp)
- (★★★☆) [文字列の相異なる連続部分文字列の個数 (O(N))](https://github.com/drken1215/algorithm/blob/master/String/num_of_substrings.cpp)


　

<a name="tr"></a>
# 木 (Tree)
木上のクエリに答えるデータ構造や、木に関する問題を解くアルゴリズムの実装です

## 木

- (★★★☆) [木の走査 (部分木サイズ, LCA, Euler tour など)](https://github.com/drken1215/algorithm/blob/master/Tree/run_tree.cpp)
- (★★★☆) [木の直径](https://github.com/drken1215/algorithm/blob/master/Tree/diameter.cpp)
- (★★★☆) 木の重心
- (★★★★) 木の Distance Frequency Table

## 木 DP

- (★★☆☆) 木 DP
- (★★★☆) 二乗の木 DP (俗称)
- (★★★☆) [全方位木 DP](https://github.com/drken1215/algorithm/blob/master/Tree/rerooting.cpp)
- (★★★☆) [全方位木 DP (辺に重みがある場合)](https://github.com/drken1215/algorithm/blob/master/Tree/rerooting_with_edge.cpp)

## LCA

- (★★★☆) [LCA (ダブリング)](https://github.com/drken1215/algorithm/blob/master/Tree/lca_by_doubling.cpp)
- (★★★☆) [LCA (Euler Tour)](https://github.com/drken1215/algorithm/blob/master/Tree/lca_euler_tour.cpp)
- (★★★☆) [LCA (HL 分解)](https://github.com/drken1215/algorithm/blob/master/Tree/lca_heavy_light_decomposition.cpp)

## 木上のクエリ処理

- (★★★☆) [Euler Tour](https://github.com/drken1215/algorithm/blob/master/Tree/euler_tour.cpp)
- (★★★☆) [HL 分解](https://github.com/drken1215/algorithm/blob/master/Tree/heavy_light_decomposition.cpp)
- (★★★☆) [重心分解](https://github.com/drken1215/algorithm/blob/master/Tree/tree_centroid_decomposition.cpp)
- (★★★☆) DSU on Tree
- (★★★★) Link-Cut 木
- (★★★★) toptree

## その他の問題

- (★★★☆) [強平衡二分木の Distance Frequency Table](https://github.com/drken1215/algorithm/blob/master/Tree/find_various_values_of_binary_tree.cpp)
- (★★★★) Level Ancester


　
<a name="ot"></a>
# その他 (OTHERS)
その他のアルゴリズムです

## 入出力

- (★★★★) [Fast IO](https://github.com/drken1215/algorithm/blob/master/Others/fast_io.cpp)

## グリッド

- (★☆☆☆) [グリッドの 4 近傍, 8 近傍](https://github.com/drken1215/algorithm/blob/master/Others/grid_neighbors.cpp)
- (★☆☆☆) [ハニカムの 6 近傍](https://github.com/drken1215/algorithm/blob/master/Others/honeycomb_neighbors.cpp)

## その他

- (★★★★) [ACL](https://github.com/drken1215/algorithm/blob/master/Others/acl.cpp)
- (★★☆☆) [XorShift, ランダムシャッフル](https://github.com/drken1215/algorithm/blob/master/Others/xorshift.cpp)
- (★★☆☆) [タイマー](https://github.com/drken1215/algorithm/blob/master/Others/timer.cpp)
- (★★☆☆) [サイコロ](https://github.com/drken1215/algorithm/blob/master/Others/dice.cpp)
- (★★☆☆) [曜日](https://github.com/drken1215/algorithm/blob/master/Others/day_of_the_week.cpp)
- (★★★☆) 四面体 (AOJ 2060)



　

# コードテンプレート

- [AtCoder 用のコードテンプレート](https://github.com/drken1215/algorithm/blob/master/template_atcoder.cpp)
- [Codeforces 用のコードテンプレート（65536 字以内)](https://github.com/drken1215/algorithm/blob/master/template_codeforces.cpp)
- [AOJ 用のコードテンプレート（30000 字以内)](https://github.com/drken1215/algorithm/blob/master/template_aoj.cpp)
- [最小限のコードテンプレート](https://github.com/drken1215/algorithm/blob/master/template_minimum.cpp)


　

# License
These codes are licensed under CC0.
[![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png "CC0")](http://creativecommons.org/publicdomain/zero/1.0/deed.ja)
