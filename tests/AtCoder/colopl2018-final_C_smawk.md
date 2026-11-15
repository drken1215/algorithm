---
title: COLOCON -Colopl programming contest 2018- Final C - スペースエクスプローラー高橋君 (NoviSteps 3D)
documentation_of: colopl2018-final_C_smawk.cpp
---

## 問題概要
長さ $N$ の数列 $a_1, a_2, \dots, a_N$ が与えられる。
各 $j = 1, 2, \dots, N$ に対して、

$\min_{1 \le i \le N} ((j - i)^2 + a_{i})$

を求めよ。


## 制約

- $1 \le N \le 2 \times 10^{5}$


## 解法メモ
$N \times N$ 行列 $A$ を

$A_{j, i} = (x - i)^2 + a_{i}$

によって定めると、行列 $A$ の各行についての最小値を求める問題となる。

$A$ は Monge である。よって、Monotone Minima 法 (計算量 $O(N \log N)$) や SMAWK 法 (計算量 $O(N)$) で解くことができる。
