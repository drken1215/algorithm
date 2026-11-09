---
title: Monge 性 & 単調性を満たす区間 DP の Knuth-Yao Speedup (in O(N^2))
documentation_of: ./knuth_yao.hpp
---

## 概要
関数 $f$ が Monge であり、かつ区間単調性をもつとき、

`dp[i][j] = min_{i <= k < j} (dp[i][k] + dp[k][j]) + f(i, j))`

という区間 DP を、$O(N^2)$ の計算量で実行する。

（この区間 DP は、一般的には $O(N^3)$ の計算量を要する形である）


## 問題例

- [KUPC 2012 J - 刺身](https://atcoder.jp/contests/kupc2012/tasks/kupc2012_10)

## 参考文献

- tatyam: [Monge の手引き書](https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu)
- 週刊 spaghetti_source: [Knuth-Yao speedup](https://topcoder-g-hatena-ne-jp.jag-icpc.org/spaghetti_source/20120915/1347668163.html)
