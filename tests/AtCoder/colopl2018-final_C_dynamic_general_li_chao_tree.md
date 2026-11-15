---
title: COLOCON -Colopl programming contest 2018- Final C - スペースエクスプローラー高橋君 (NoviSteps 3D)
documentation_of: colopl2018-final_C_dynamic_general_li_chao_tree.cpp
---

## 問題概要
長さ $N$ の数列 $a_1, a_2, \dots, a_N$ が与えられる。
各 $j = 1, 2, \dots, N$ に対して、

$\min_{1 \le i \le N} ((j - i)^2 + a_{i})$

を求めよ。


## 制約

- $1 \le N \le 2 \times 10^{5}$


## 解法メモ
関数群 $f_1, f_2, \dots, f_N$ を

$f_{i}(x) = (x - i)^2 + a_{i}$

によって定める。このとき、各 $j = 1, 2, \dots, N$ に対して、

$\min_{1 \le i \le N} f_{i}(j)$

を求める問題となる。

関数群 $f_1, f_2, \dots, f_N$ は Monge である。よって、Li Chao Tree を用いると、

- 関数 $f_1, f_2, \dots, f_N$ を Li Chao Tree に挿入する
- $\min_{1 \le i \le N} f_{i}(j)$ を求める

によって解くことができる。計算量は $O(N \log N)$ である。
