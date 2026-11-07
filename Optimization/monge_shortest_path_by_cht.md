---
title: Monge 単一始点最短路問題 (by Convex Hull Trick, in O(N log N))
documentation_of: ./monge_shortest_path_by_cht.hpp
---

## 概要
Convex Hull Trick (General Li Chao Tree) を用いた単一始点最短路アルゴリズム。

DP が

`dp[j] = min_{i < j} (dp[i] + func(i, j))`

の形で表されるとする。ここで、

`g(j, i) = dp[i] + func(i, j)` 

が PSC (Pairwise Single Crossing) 性を満たすとき、Convex Hull Trick (General Li Chao Tree) を用いて、
$O(N \log N)$ の計算量で DP テーブルを求めることができる。

とくに、`func(l, r)` が $0 \le l, r \le N$ を満たすすべての $l, r$ に対して Monge 性、または anti-Monge 性を満たすとき、上記の仮定は満たされる。

## 仕様

- `template<class VAL, class FUNC> vector<pair<VAL, int>> MongeShortestPathByCHT(int N, const FUNC &f)`
  - コスト関数が $f$ であるような、次のグラフの単一始点最短路を求める
    - 頂点: $0, 1, 2, \dots, N$
    - 辺: $(l, r)$ ($0 \le l < r \le N$)
- `template<class VAL, class FUNC> vector<pair<VAL, int>> MongeShortestPathByCHT(int N, const FUNC &f, const vector<int> &start)`
  - コスト関数が $f$ であるような、次のグラフの単一始点最短路を求める
    - 頂点: $0, 1, 2, \dots, N$
    - 辺: $(l, r)$ ($0 \le l < r \le N$, $r \ge$ $\mathrm{start}[l]$) ($\mathrm{start}[l] \ge l + 1$ でなければならない)

## 計算量
$O(N \log N)$

## 注意点
問題設定におけるコスト関数 $f(l, r)$ の定義域は $0 \le l < r \le N$ であるが、$l \ge r$ であるような $l, r$ に対しても定義されなければならない。

その上で、`g(j, i) = dp[i] + func(i, j)` が PSC 性を満たすことを要請する。

## 問題例

- [AtCoder ABC 480 G - Zonal Score Minimization](https://atcoder.jp/contests/abc480/tasks/abc480_g)
