---
title: Monge 単一始点 d-辺最短路の d = 1, 2, ..., D における列挙
documentation_of: ./monge_enum_d_shortest_paths.hpp
---

## 概要

- 頂点: $0, 1, 2, \dots, N$
- 辺: $(l, r)$ ($0 \le l < r \le N$), コスト $f(l, r)$

であるような DAG を考える。$f(l, r)$ が Monge であるとする。

このとき、$d = 1, 2, \dots, D$ に対して、頂点 $0$ から各頂点への最短パスのうち、パスに含まれる辺の本数がちょうど $d$ 本であるようなものを求める。


### DP で解釈
次の DP を高速化するアルゴリズムであるとも言える。

`dp[j][d] = min_{j < i} (dp[i][d - 1] + f(i, j))`

これは通常 $O(N^2 D)$ の計算量を要するが、$f$ が Monge ならば高速に実行できる。

- SMAWK 法を用いると $O(ND)$ の計算量で実行できる
- Monotone Minima 法を用いると $O(DN \log N)$ の計算量で実行できる 
  - 辺コストをスライドアクセスで求めたい場合は SMAWK 法は使えないので、これが重宝する


## 仕様

### ランダムアクセス版


### スライドアクセス版


## 問題例

- [yukicoder No.952 危険な火薬庫](https://yukicoder.me/problems/no/952)
- [Monoxer Programming Contest for Engineers K - Coupon](https://mofecoder.com/contests/monoxercon_202508/tasks/monoxercon_202508_k)
- [Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem](https://codeforces.com/contest/868/problem/F) (スライドアクセス)

## 参考文献

- tatyam: [Monge の手引き書](https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu)
- noshi: [簡易版 LARSCH Algorithm](https://noshi91.hatenablog.com/entry/2023/02/18/005856)
- 週刊 spaghetti_source: [Knuth-Yao speedup](https://topcoder-g-hatena-ne-jp.jag-icpc.org/spaghetti_source/20120915/1347668163.html)
