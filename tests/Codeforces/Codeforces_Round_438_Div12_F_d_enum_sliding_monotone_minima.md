---
title: 分割統治 Monotone Minima 法において、辺コストをスライドアクセスで計算する例
documentation_of: ./Codeforces_Round_438_Div12_F_d_enum_sliding_monotone_minima.cpp
---

## 概要
DP が

- `dp[j] = min_{i < j} (dp[i] + f(i, j))`
- `dp[j][d] = min_{j < i} (dp[i][d - 1] + f(i, j))`

といった形で表されるとする。ここで、$f$ が Monge であるならば、上記の DP は高速に実行できる。

ここで、`f(i, j)` が両端の伸縮にともなう差分更新によって求められる (スライドアクセス) 場合には、

- Monotone Minima 法は使える
- SMAWK 法は使えない
- DP 更新式のうち前者につては、noshi's 簡易 LARSCH 法は使える

ここでは、後者の DP 更新について、Monotone Minima 法による実装例を示す。

## 問題例

- [Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem](https://codeforces.com/contest/868/problem/F)

## 参考文献

- noshi: [簡易版 LARSCH Algorithm](https://noshi91.hatenablog.com/entry/2023/02/18/005856)
