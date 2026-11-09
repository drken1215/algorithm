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
template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f)
```

構造体 `DDSMAWK<VAL>` を宣言して、関数 `solve(N, f)` を呼び出す。

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
- 週刊 spaghetti_source: [Totally Monotone Matrix Searching (SMAWK algorithm)](https://topcoder-g-hatena-ne-jp.jag-icpc.org/spaghetti_source/20120923/1348327542.html)
- tmaehara: [動的計画法高速化テクニック（オフライン・オンライン変換）](https://qiita.com/tmaehara/items/0687af2cfb807cde7860)
