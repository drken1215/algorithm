---
title: Monge 単一始点最短路問題 (by Convex Hull Trick, in O(N log N))
documentation_of: ./monge_shortest_path_by_cht.hpp
---

## 概要
[General Li Chao Tree](https://drken1215.github.io/algorithm/Optimization/general_li_chao_tree.hpp.html) を用いた Monge 単一始点最短路アルゴリズム。

頂点数 $N + 1$ の DAG で、辺 $(l, r)$ ($0 \le l < r \le N$) のコスト $f(l, r)$ が Monge 性を満たすとき、Convex Hull Trick (Li Chao Tree) を用いて、頂点 0 から各頂点への最短路を求めることができる。


## 仕様

- `template<class VAL, class FUNC> vector<pair<VAL, int>> MongeShortestPathByCHT(int N, const FUNC &f)`
  - コスト関数が $f$ であるような、次のグラフの単一始点最短路を求める
    - 頂点: $0, 1, 2, \dots, N$
    - 辺: $(l, r)$ ($0 \le l < r \le N$)
- `template<class VAL, class FUNC> vector<pair<VAL, int>> MongeShortestPathByCHT(int N, const FUNC &f, const vector<int> &start)`
  - コスト関数が $f$ であるような、次のグラフの単一始点最短路を求める
    - 頂点: $0, 1, 2, \dots, N$
    - 辺: $(l, r)$ ($r \ge$ `start`[$l$]) (`start`[$l$] $\le l + 1$ でなければならない)

## 計算量
$O(N \log N)$

## 制約
問題設定におけるコスト関数 $f(l, r)$ の定義域は $0 \le l < r \le N$ であるが、$l \ge r$ を満たす $(l, r)$ に対しても定義され、全体で Monge 性を満たす必要がある

## 問題例

- [AtCoder ABC 480 G - Zonal Score Minimization](https://atcoder.jp/contests/abc480/tasks/abc480_g)
