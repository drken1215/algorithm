---
title: TM 単一始点最短路問題 (D&D SMAWK 法, in O(N log N))
documentation_of: ./dd_smawk.hpp
---

## 概要
SMAWK 法に対して、オンライン・オフライン変換を施して得られる分割統治 SMAWK 法。単一始点最短路問題を解く。

DP が

`dp[j] = min_{i < j} (dp[i] + func(i, j))`

の形で表されるとする。ここで、

`g(j, i) = dp[i] + func(i, j)` 

が Totally Monotone であるならば、上記の DP が分割統治 SMAWK 法によって、$O(N \log N)$ の計算量で DP テーブルを求めることができる。

## 注意点

- `func(i, j)` が Totally Monotone であっても、`g(j, i) = dp[i] + func(i, j)` が Totally Monotone であるとは限らない
- `func(i, j)` が Monge であるとき、`g(j, i) = dp[i] + func(i, j)` も Monge である

上記のことから、多くの場合、コストが Monge 性を満たすときに適用できるアルゴリズムであるとされる。

## 仕様

```cpp
template<class VAL, class FUNC> vector<pair<VAL, int>> DDSMAWK(int N, const FUNC &f)
```

- 整数 `N`: 頂点番号が $0, 1, 2, \dots, N$ であるとする (頂点数は $N + 1$)
- 関数 `f(l, r)`: 辺 $(l, r)$ のコストを返す関数
- 返り値: 受け取る変数を `dp` として
  - `dp[v].first`: 頂点 $0$ から頂点 $v$ への最短路の重み
  - `dp[v].second`: 頂点 $0$ から頂点 $v$ への最短路において、$v$ に至る直前の頂点の番号 (存在しない場合は -1)


## 問題例

- [AtCoder EDPC Z - Frog 3](https://atcoder.jp/contests/dp/tasks/dp_z)
- [Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry](https://codeforces.com/contest/319/problem/C)
- [yukicoder No.705 ゴミ拾い Hard](https://yukicoder.me/problems/no/705)

## 参考文献

- tatyam: [Monge の手引き書](https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu)
