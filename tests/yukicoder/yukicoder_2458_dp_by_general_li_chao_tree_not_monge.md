---
title: yukicoder No.2458 Line Up Charged Balls
documentation_of: ./yukicoder_2458_dp_by_general_li_chao_tree_not_monge.test.cpp
---

## 問題概要
長さ $N$ の数列 $Q_1, \dots, Q_N$ が与えられる。

この数列から順序を変えずに 2 個以上を抜き取って得られる数列を $q_1, \dots, q_K$ とするとき、

$q_1 q_2 + q_2 q_3 + \dots + q_{K-1} q_{K}$

の最大値を求めよ。


### 制約

- $3 \le N \le 3 \times 10^{5}$
- $-10^{5} \le Q_{i} \le 10^{5}$


## 解法メモ
$f(i, j) = Q_{i} Q_{j}$ 

として、次の DP で最適解が求められる問題である。

`dp[j] = min_{j < i} (dp[i] + f(i, j))`

今回は 

- $f(i, j)$ は Monge ではない
- $g_{i}(j) = Q_{i}Q_{j} + \mathrm{dp[i]}$ とおいたとき、関数群 $g_1, g_2, \dots$ は Pairwise Single Crossing でもない
- しかし、$g_{i}(x) = Q_{i}x + \mathrm{dp[i]}$ とおいたとき、関数群 $g_1, g_2, \dots$ は Pairwise Single Crossing である

よって、Convex Hull Trick (Li Chao Tree) によって解くことができる。
