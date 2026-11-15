---
title: CCOLOCON 2018 Final C - スペースエクスプローラー高橋君 (単調 Convex Hull Trick を用いる方法) (NoviSteps 3D)
documentation_of: colopl2018-final_C_monotone_convex_hull_trick.cpp
---

## 問題概要
長さ $N$ の数列 $a_1, a_2, \dots, a_N$ が与えられる。
各 $j = 1, 2, \dots, N$ に対して、

$\min_{1 \le i \le N} ((j - i)^2 + a_{i})$

を求めよ。


## 制約

- $1 \le N \le 2 \times 10^{5}$


## 解法メモ
$((j - i)^2 + a_{i})$ = (-2ix + (i^2 + a_{i})) + x^2$

であるから、直線群 $f_1, f_2, \dots, f_N$ を

$f_{i}(x) = -2i x + (i^2 + a_{i})$

によって定めると、各 $j = 1, 2, \dots, N$ に対して、

$j^2 + \min_{1 \le i \le N} f_{i}(j)$

を求める問題となる。ここで、傾き $-2i$ が単調減少であることから、傾きが単調減少であることを仮定した Deque ベースの Convex Hull Trick によって、

- 関数 $f_1, f_2, \dots, f_N$ をこの順に挿入する (挿入の計算量は $O(1)$)
- $\min_{1 \le i \le N} f_{i}(j)$ を求める ($j$ も単調増加であるため、クエリ取得の計算量も $O(1)$)

によって解くことができる。計算量は $O(N)$ である。
