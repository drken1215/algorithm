---
title: Monge 行列の min-plus 合成 (in O(N^2))
documentation_of: ./monge_min_plus_composition.hpp
---

## 概要
上三角 $(N+1) \times (N+1)$ Monge 行列 $A, B$ に対して

`C[i][j] = min_{i <= k < j} (A[i][k] + A[k][j])`

を高速に実行する。

いわゆる「斜めに埋めるテク」を用いることで、$O(N^2)$ の計算量で実行できる。

## 問題例

## 参考文献

- tatyam: [Monge の手引き書](https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu)
- 週刊 spaghetti_source: [Knuth-Yao speedup](https://topcoder-g-hatena-ne-jp.jag-icpc.org/spaghetti_source/20120915/1347668163.html)
