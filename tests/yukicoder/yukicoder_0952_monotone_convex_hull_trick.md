---
title: yukicoder No.952 危険な火薬庫 (NoviSteps 3D+)
documentation_of: ./yukicoder_0952_monotone_convex_hull_trick.test.cpp
---

## 問題概要
数列 $A_1, A_2, \dots, A_N$ が与えられる。$k = 1, 2, \dots, N$ に対して、次の問いに答えよ。

- 数列を $k$ 個の区間に分割したい
- 各区間のスコアは、各区間に属する数値の総和の 2 乗で与えられる

各区間のスコアの総和の最大値を求めよ。



## 制約

- $1 \le N \le 3000$


## 解法メモ
数列 $A$ の累積和を $S$ とすると、

+ $f(i, j) = (S_{j} - S_{i})^{2}$

として、次の DP が立てられる。

`dp[j][k] = min_{j < i} (dp[i][k-1] + f(i, j))`

(`dp[i][k]` := 先頭から $i$ 個のみを考えて、それを $k$ 個の区間に分割するときのスコアの最大値)

$f(i, j)$ は Monge であるから、次のような解法によって解くことができる。

- Monge 単一始点 $d$-辺最短路の $d = 1, 2, \dots, D$ における列挙 (計算量 $O(N^2)$)
- General Li Chao Tree を用いる (計算量 $O(N^2 \log N)$)
- $f$ を展開して、直線の Convex Hull Trick に帰着して解く (傾きの単調性を活かすことで、計算量は $O(N^2)$ となる)

本コードは、「$f$ を展開して、直線の Convex Hull Trick に帰着して解く」ものである。
