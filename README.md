# 様々なアルゴリズム
データ構造や数論的アルゴリズムまで、様々な分野のアルゴリズムたちを C++14 で実装しています。  
競技プログラミングなどの場面においてライブラリとして活用することを念頭に置いています。

# DataStructure
各種データ構造の実装です

- union-find 木
  - [union-find 木 (基本)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_tree_simple.cpp)
  - [union-find 木 (rank つき)](https://github.com/drken1215/algorithm/blob/master/DataStructure/union_find_tree.cpp)
  - [重みつき union-find 木](https://github.com/drken1215/algorithm/blob/master/DataStructure/weighted_union_find_tree.cpp)
- BIT
  - [BIT (基本)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree.cpp)
  - [BIT 上二分探索 (k 番目の要素を求める)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_search_on_BIT.cpp)
  - [BIT (区間加算, 区間和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_RAQ.cpp)
  - [二次元 BIT](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_2D.cpp)
  - [二次元 BIT (領域加算, 領域和取得に両対応)](https://github.com/drken1215/algorithm/blob/master/DataStructure/binary_indexed_tree_2D_RAQ.cpp)
- RMQ
  - [RMQ (セグメントツリー)](https://github.com/drken1215/algorithm/blob/master/DataStructure/range_minimum_query.cpp)
  - [RMQ (sparse table)](https://github.com/drken1215/algorithm/blob/master/DataStructure/sparse_table.cpp)
- セグメントツリー
  - [セグメントツリーテンプレ](https://github.com/drken1215/algorithm/blob/master/DataStructure/segment_tree.cpp)
- 平衡二分探索木
  - [RBST](https://github.com/drken1215/algorithm/blob/master/DataStructure/randomized_binary_search_tree.cpp)


# DataStructureOnTree
ツリー上のクエリ処理のためのデータ構造たちの実装です

- [HL 分解](https://github.com/drken1215/algorithm/blob/master/DataStructureOnTree/heavy_light_decomposition.cpp)


# DP
定型的な動的計画法やその他の処理です

- [いもす法 (俗称)](https://github.com/drken1215/algorithm/blob/master/DP/imos.cpp)
- [転倒数](https://github.com/drken1215/algorithm/blob/master/DP/inversion_number.cpp)


# Geometry
幾何ライブラリです

- [全部乗せ](https://github.com/drken1215/algorithm/blob/master/Geometry/All.cpp)
- [基本要素 (点, 線分, 円)](https://github.com/drken1215/algorithm/blob/master/Geometry/BasicElements.cpp)
- [射影](https://github.com/drken1215/algorithm/blob/master/Geometry/Projection.cpp)


# GraphNetworkFlow
グラフネットワークフロー関連のアルゴリズムです

- [最大流 (Dinic 法)](https://github.com/drken1215/algorithm/blob/master/GraphNetworkFlow/max_flow_dinic.cpp)


# GraphTheory
グラフ理論全般のアルゴリズムです

- [最大安定集合問題 (O(1.381^n))](https://github.com/drken1215/algorithm/blob/master/GraphTheory/maximum_stable_set.cpp)


# MathAlgebra
行列計算など代数的計算に関するアルゴリズムです

- [行列]()
- [行列累乗]()
- [連立一次方程式 (実数)]()
- [連立一次方程式 (mod. p)]()
- [連立一次方程式 (bitset)]()


# MathNumberTheory
整数論的アルゴリズムたちです、mod 演算も含みます

- [mod 演算 (逆元、累乗、二項係数)](https://github.com/drken1215/algorithm/blob/master/MathNumberTheory/mod.cpp)


# String
文字列アルゴリズムです

- [suffix array](https://github.com/drken1215/algorithm/blob/master/String/suffix_array.cpp)


# Others
その他のアルゴリズムです



