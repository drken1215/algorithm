---
title: General Li Chao Tree (一般関数対応)
documentation_of: ./general_li_chao_tree.hpp
---

## 概要
一般関数 (Pairwise Single Crossing) に対応した Li Chao Tree

### Dynamic 版 (座標圧縮しない)

- コンストラクタ: `DynamicGeneralLiChaoTree(X xl, X xr, VAL inf = numeric_limits<VAL>::max() / 2)`
  - $x$ の範囲が `[xl, xr)` (この幅を $X$ とする)
- 関数の挿入
  - `insert(f)`: 関数 $f$ を挿入する。計算量 $O(\log X)$ 
  - `insert_segment(f, l, r)`: 関数 $f$ を区間 $[l, r)$ の範囲内で挿入する。計算量 $O(\log^{2} X)$
- 1 点取得
  - `query(x)`: $\min_{i} f_{i}(x)$ を返す。計算量 $O(\log X)$

### Static 版 (座標圧縮するか、定義域の範囲が小さい)



## 制約

- 一般関数 $f_1, f_2, \dots, f_N$ が Pairwise Single Crossing であること

## 問題例

- [Library Checker - Line Add Get Min](https://judge.yosupo.jp/problem/line_add_get_min) (直線を管理)
- [Library Checker - Segment Add Get Min](https://judge.yosupo.jp/problem/segment_add_get_min) (直線を管理, 区間挿入あり)
- [COLOCON 2018 Final C - スペースエクスプローラー高橋君](https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c) (一般関数)
- [AtCoder EDPC Z - Frog 3](https://atcoder.jp/contests/dp/tasks/dp_z) (Monge 単一始点最短路)
- [Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry](https://codeforces.com/contest/319/problem/C) (Monge 単一始点最短路)
- [yukicoder No.705 ゴミ拾い Hard](https://yukicoder.me/problems/no/705) (Monge 単一始点最短路)
- [AtCoder ABC 480 G - Zonal Score Minimization](https://atcoder.jp/contests/abc480/tasks/abc480_g) (anti-Monge 単一始点最短路)

## 参考文献

- Algorithms for Competitive Programming: [Convex hull trick and Li Chao tree](https://cp-algorithms.com/geometry/convex_hull_trick.html)
- robert1003: [A Simple Introduction to Li-Chao Segment Tree](https://robert1003.github.io/2020/02/06/li-chao-segment-tree.html)
- rama_pang: [[Tutorial] Li Chao Tree Extended](https://codeforces.com/blog/entry/86731)
