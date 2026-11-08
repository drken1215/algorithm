---
title: Monge 単一始点最短路問題 (by LARSCH 法, in O(N))
documentation_of: ./larsch.hpp
---

## 概要
DP が

`dp[j] = min_{i < j} (dp[i] + func(i, j))`

の形で表されるとする。ここで、`func(i, j)` が Monge であるとき、

`g(j, i) = dp[i] + func(i, j)` 

も Monge であり、上記の DP は LARSCH 法を用いて $O(N)$ の計算量で実行できる。


## 仕様

```cpp
template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f)
```

構造体 `LARSCH<VAL>` を宣言して、関数 `solve(N, f)` を呼び出す。

- 整数 `N`: 頂点番号が $0, 1, 2, \dots, N$ であるとする (頂点数は $N + 1$)
- 関数 `f(l, r)`: 辺 $(l, r)$ のコストを返す関数
- 返り値: 受け取る変数を `res` として
  - `res[v].first`: 頂点 $0$ から頂点 $v$ への最短路の重み
  - `res[v].second`: 頂点 $0$ から頂点 $v$ への最短路において、$v$ に至る直前の頂点の番号 (存在しない場合は -1)
- また、構造体のメンバ変数にも結果が格納される
  - `dp[v]`: 頂点 $0$ から頂点 $v$ への最短路の重み
  - `cnt[v]`: 頂点 $0$ から頂点 $v$ への最短路に含まれる辺の本数
  - `prev[v]`: 頂点 $0$ から頂点 $v$ への最短路において、$v$ に至る直前の頂点の番号 (存在しない場合は -1)


```cpp
vector<int> reconstruct()
```

- `solve(N, f)` 実行後に呼び出すことで、頂点 $0$ から頂点 $N$ への最短路が復元される


## 問題例

- [AtCoder EDPC Z - Frog 3](https://atcoder.jp/contests/dp/tasks/dp_z)
- [Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry](https://codeforces.com/contest/319/problem/C)
- [yukicoder No.705 ゴミ拾い Hard](https://yukicoder.me/problems/no/705)

## 参考文献

- tatyam: [Monge の手引き書](https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu)
