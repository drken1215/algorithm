# 様々なアルゴリズム
データ構造や数論的アルゴリズムまで、様々な分野のアルゴリズムたちを C++14 で実装しています。  
プログラミングコンテストなどにおいて「実装例」または「ライブラリ」として使用することを念頭に置いています。

  

# DataStructure
各種データ構造の実装です

#### Union-Find 木

- [Union-Find 木 (union by size)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_tree_size.cpp)
- [Union-Find 木 (union by rank)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_tree_rank.cpp)
- [重みつき Union-Find 木](https://github.com/drken1215/algorithm/blob/master/DataStructure/weighted_union_find_tree.cpp)
- [部分永続 Union-Find 木](https://github.com/drken1215/algorithm/blob/master/DataStructure/partially_persistent_union_find_tree.cpp)

#### セグメント木

- [セグメント木](https://github.com/drken1215/algorithm/blob/master/DataStructure/segment_tree.cpp)
- [セグメント木 (遅延評価)](https://github.com/drken1215/algorithm/blob/master/DataStructure/segment_tree_delay.cpp)
- [Starry Sky 木 (俗称)](https://github.com/drken1215/algorithm/blob/master/DataStructure/starry_sky_tree.cpp)
- 二次元セグメント木

#### BIT

- [BIT](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree.cpp)
- [BIT 上二分探索 (k 番目の要素を求める)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_search_on_BIT.cpp)
- [BIT (区間加算, 区間和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_RAQ.cpp)
- [二次元 BIT](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_2D.cpp)
- [二次元 BIT (領域加算, 領域和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_2D_RAQ.cpp)

#### RMQ

- [RMQ (セグメント木)](https://github.com/drken1215/algorithm/blob/master/DataStructure/range_minimum_query.cpp)
- [RMQ (Sparse Table)](https://github.com/drken1215/algorithm/blob/master/DataStructure/sparse_table.cpp)

#### 平衡二分探索木

- [RBST](https://github.com/drken1215/algorithm/blob/master/DataStructure/randomized_binary_search_tree.cpp)
- Splay 木
- 赤黒木

#### 永続データ構造

- 永続配列
- 完全永続 Union-Find 木
- 永続セグメント木
- 永続赤黒木

#### その他

- [Disjoint Sparse Table](https://github.com/drken1215/algorithm/blob/master/DataStructure/disjoint_sparse_table.cpp)
- 並列二分探索
- Skew Heap
- Wavelet 木




  
# DataStructureOnTree
ツリー上のクエリ処理のためのデータ構造たちの実装です

#### LCA

- [LCA (ダブリング)](https://github.com/drken1215/algorithm/blob/master/DataStructureOnTree/LCA_doubling.cpp)
- LCA (Euler Tour)
- LCA (HL 分解)

#### テクニック

- Euler Tour
- 重心分解
- [HL 分解](https://github.com/drken1215/algorithm/blob/master/DataStructureOnTree/heavy_light_decomposition.cpp)
- Link-Cut 木




  
# DP
定型的な動的計画法やその他の処理です

#### いもす法 (俗称)

- [いもす法 (俗称)](https://github.com/drken1215/algorithm/blob/master/DP/imos.cpp)
- [二次元いもす法 (俗称)](https://github.com/drken1215/algorithm/blob/master/DP/imos_2D.cpp)
- 三次元いもす法 (俗称)

#### 典型的 DP

- [転倒数](https://github.com/drken1215/algorithm/blob/master/DP/inversion_number.cpp)
- LIS
- LCS
- 編集距離
- ヒストグラム長方形面積最大化
- 最適二分探索木

#### DP 高速化テクニック

- 累積和
- スライド最小値
- Convex Hull Trick
- Monotone Minima
- Monge




  
# Geometry
幾何ライブラリです

- [全部乗せ](https://github.com/drken1215/algorithm/blob/master/Geometry/All.cpp)
- [基本要素 (点, 線分, 円)](https://github.com/drken1215/algorithm/blob/master/Geometry/BasicElements.cpp)

#### 射影, 距離, 交差判定

- [射影](https://github.com/drken1215/algorithm/blob/master/Geometry/Projection.cpp)
- 距離
- 交点
- 線分アレンジメント

#### 多角形

- 多角形の面積
- 点と多角形の包含判定
- 凸包
- 凸多角形の切断
- 凸多角形の直径
- ボロノイ図 (単純ver, O(n^2))

#### 接線

- 接線

#### 三次元幾何

- 基本要素

#### その他

- 最近点対
- 最近円対
- 双対変換
- kd 木



  
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

#### マッチング

- 二部マッチング (Hopcroft-Karp 法)
- 重みつき二部マッチング (Hungarian 法)
- 一般グラフの最大マッチング (Edmonds 法)
- 一般グラフの最大マッチング (行列補間)




  
# GraphTheory
グラフ理論全般のアルゴリズムです

#### 連結成分分解

- 連結成分分解
- 強連結成分分解
- 二重辺連結成分分解
- 二重点連結成分分解

#### ツリー

- ツリーの直径
- ツリーの重心

#### 最短路

- 単一始点最短路 (Dijkstra 法, 正辺のみ)
- 単一始点最短路 (Bellman-Ford 法, 負辺対応)
- 全頂点間最短路 (Floyd-Warshall 法)
- k-最短路
- SPFA

#### その他

- 最小全域木問題 (Kruskal 法)
- Euler 路
- [最大安定集合問題 (O(1.381^n))](https://github.com/drken1215/algorithm/blob/master/GraphTheory/maximum_stable_set.cpp)




  
# MathAlgebra
行列計算など代数的計算に関するアルゴリズムです

#### 行列

- [行列 (基本演算)](https://github.com/drken1215/algorithm/blob/master/MathAlgebra/matrix.cpp)
- 行列累乗 (実数)
- 行列累乗 (mod. p)
- 行列累乗 (binary)
- 連立一次方程式 (実数)
- 連立一次方程式 (mod. p)
- 連立一次方程式 (binary)
- Toeplitz 行列

#### 多項式, 方程式

- 二次方程式
- 多項式 (実数係数)
- 多項式 (mod. p 係数)
- きたまさ法 (俗称)
- きたまさ法 with FFT (俗称)

#### 畳み込み計算

- FFT (高速フーリエ変換)
- NTT (高速剰余変換)
- 高速アダマール変換
- 高速ゼータ変換
- 高速メビウス変換
- Karatsuba 法

#### 最適化

- 三分探索法
- Newton 法
- 単体法




  
# MathCombinatorics
組合せ論的アルゴリズムたちです

#### mod, 二項係数

- [mod 演算 (全部乗せ)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod.cpp)
- [mod 演算 (累乗)](https://github.com/drken1215/algorithm/blob/master/MathCombinatorics/mod_power.cpp)
                  

#### 様々な数

- カタラン数
- 分割数
- スターリング数
- ベル
#### ソート

- Radix ソート

#### マトロイド

- マトロイド上の Greedy 法
- マトロイド交差




  
# MathNumberTheory
整数論的アルゴリズムたちです

#### 約数, 倍数

- [最大公約数 (Euclid の互除法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/GCD.cpp)
- [最小公倍数](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/LCM.cpp)
- [拡張 Euclid の互除法](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/extended_GCD.cpp)

#### 素数

- [素数判定](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/is_prime.cpp)
- [素因数分解](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/prime_factorization.cpp)
- [確率的素数判定 (Miller-Rabin 法)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/is_prime_Miller_Rabin.cpp)
- 確率的素因数分解
- エラトステネスの篩
- エラトステネスの区間篩
- アトキンの篩

#### 方程式

- 中国剰余定理
- 中国剰余定理 (Garner 法)
- 連立一次合同方程式
- ペル方程式
- 離散対数
- 平方剰余

#### 有理数

- 有理数
- Stern-Brocot 木

#### その他
  
- 多倍長整数




  
# String
文字列アルゴリズムです

#### 文字列検索

- ローリングハッシュ
- 単一パターン検索 (KMP 法)
- 単一パターン検索 (Boyer-Moore 法)
- 複数パターン検索 (Aho-Corasick 法)

#### 構文解析

- LL(1) 再帰降下パーサ

#### 文字列系データ構造

- Trie 木
- [Suffix Array](https://github.com/drken1215/algorithm/blob/master/String/suffix_array.cpp)
- Z 法
- Mo 法
- Manacher 法




  
# Others
その他のアルゴリズムです

- サイコロ
- 曜日
- α-β 探索
- A*
- IDA*
- 2-SAT










