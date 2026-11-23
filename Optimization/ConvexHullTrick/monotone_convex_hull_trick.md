---
title: Convex Hull Trick (直線：傾き単調, by Deque, in $O(\log N)$) (クエリも単調なら $O(1)$)
documentation_of: ./monotone_convex_hull_trick.hpp
---

## 概要
挿入する直線の傾きが**単調減少**であることを前提とした Convex Hull Trick。次のクエリに答えるデータ構造である。

----

- 【挿入】：直線 $f_{i}(x) = a_{i} x + b_{i}$ を順に挿入する (計算量はならし $O(1)$)
- 【取得】：$x$ 座標が与えられて、$\min_{i}f_{i}(x)$ の値を出力する (計算量は後述)

----

ここで、

- $a_{0} \ge a_{1} \ge a_{2} \ge \dots$

であることを制約として課す。取得クエリの $x$ 座標については単調性を課さないが、課すと計算量が改善する。

- 取得クエリの $x$ 座標が単調増加である場合：$O(1)$
- 取得クエリの $x$ 座標に単調性を課さない場合：$O(\log N)$


## 問題例

- [AtCoder EDPC Z - Frog 3](https://atcoder.jp/contests/dp/tasks/dp_z) (Monge 単一始点最短路) (NoviSteps 3D+)

## 参考文献

- えび: [Convex Hull Trick](https://hcpc-hokudai.github.io/archive/algorithm_convex_hull_trick_001.pdf)
