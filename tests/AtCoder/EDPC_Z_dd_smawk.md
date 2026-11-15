---
title: EDPC Z - Frog 3 (分割統治 SMAWK 解法) (NoviSteps 3D+)
documentation_of: EDPC_Z_dd_smawk.cpp
---

## 問題概要
$N$ 個の足場 $1, 2, \dots, N$ がある。最初、足場 $1$ にカエルがいて、以下の動きを繰り返す。

- 足場 $i$ にいるとき、足場 $i + 1, i + 2, \dots, N$ のいずれかへジャンプする
  - このとき、ジャンプ先の足場を $j$ とすると、コスト $(h_{j} - h_{i})^2 + C$ を支払う

カエルが足場 $N$ に辿り着くまでに支払うコストの最小値を求めよ。


## 制約

- $2 \le N \le 2 \times 10^{5}$
- $1 \le h_{1} < h_{2} < \dots < h_{N} \le 10^{6}$


## 解法メモ
次のような DP で解くことができる。

+ `dp[i]` := 足場 $i$ に辿り着くまでのコスト

として、

`dp[j] = min_{i < j} (dp[i] + f(i, j))`

ここで、$f(i, j) = (h_{j} - h_{i})^2 + C$ とする。

$f$ は Monge であるから、さまざまな解法で解くことができる。

+ 分割統治 Monotone Minima 法 (計算量 $O(N \log^{2} N))
+ 分割統治 SMAWK 法 (計算量 $O(N \log N))
+ noshi's 簡易 LARSCH 法 (計算量 $O(N \log N))
+ LARSCH 法 (計算量 $O(N))
+ Convex Hull Trick を用いる方法 ($h$ の単調性を用いる) (計算量 $O(N))
+ General Li Chao Tree を用いる方法 ($h$ の単調性がなくても解ける方法) (計算量 $O(N \log N))

本ページのコードは、分割統治 SMAWK 法を用いている。
