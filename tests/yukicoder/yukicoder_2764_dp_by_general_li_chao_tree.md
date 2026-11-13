---
title: yukicoder No.2764 Warp Drive Spacecraft
documentation_of: ./yukicoder_2764_dp_by_general_li_chao_tree.test.cpp
---

## 問題概要
$N$ 個の惑星 $1, 2, \dots, N$ からなる銀河がある。銀河で惑星間を移動する方法は以下の 2 通りがある。

- 航路に沿った移動 ($M$ 本の航路があり、航路 $i$ は惑星 $U_i, V_i$ 間を双方子に繋いでいて、時間 $T_i$ で移動できる)
- ワープ航法による移動 (惑星 $i, j$ 間を、時間 $W_i W_j$ で移動できる)

ワープ航法による移動は高々 1 回しか使えない。

惑星 $1$ から惑星 $N$ へとたどり着くまでの最短時間を求めよ。

### 制約

- $N, M \le 2 \times 10^{5}$


## 解法メモ
次の 2 テーブルを管理する方針をとった。

- `dp[v]`: ワープ航法を用いずに惑星 $v$ へと到達するまでの最短時間
- `dp2[v]`: ワープ航法を 1 回以下だけ用いて惑星 $v$ へと到達するまでの最短時間

次の順序で遷移を行った。

- `dp` 間の遷移: Dijkstra 法
- `dp` → `dp2` 間の遷移: Convex Hull Trick を用いる (後述)
- `dp2` 間の遷移: Dijkstra 法

`dp` → `dp2` 間の遷移については、次のように書ける。

`dp2[j] = min_{1 <= i <= N} (dp[i] + W_i W_j)`

この遷移は、次の機能を持たせた Convex Hull Trick を用いて実現できる。

- 各 $i$ に対して、$f_{i}(x) = \mathrm{dp}[i] + W_{i} x$ を挿入する
- 各 $j$ に対して、$\min_{i} f_{i}(W_j)$ を求める

Li Chao Tree が使える形である。


## 解法選択の余地
$W$ が小さい順に惑星を並び直すと、「`dp` → `dp2` 間の遷移」に登場するコスト関数は 

- 挿入する直線の傾きも、クエリも、単調性をもつような Convex Hull Trick が使える形になる
- **anti-Monge** 性を満たすため、SMAWK 法や Monotone Mimina 法を使って解くこともできる


## 本ページの解法
本ページでは、$W$ を小さい順に並び替えることはせず、Li Chao Tree を用いている。計算量は

- Dijkstra パート: $O(M \log N)$
- Convex Hull Trick パート: $O(N \log N)$

である。
