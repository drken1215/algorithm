---
title: Monge 単一始点最短路問題 (by noshi's 簡易 LARSCH 法, in O(N log N))
documentation_of: ./simple_larsch.hpp
---

## 概要
DP が

`dp[j] = min_{i < j} (dp[i] + func(i, j))`

の形で表されるとする。ここで、`func(i, j)` が Monge であるとき、

`g(j, i) = dp[i] + func(i, j)` 

も Monge であり、上記の DP は、noshi's 簡易 LARSCH 法を用いて $O(N \log N)$ の計算量で実行できる。


### LARSCH 法との比較
上記の DP は、`func(i, j)`が Monge であるとき、LARSCH 法によって $O(N)$ の計算量でも実行できることが知られているが、これは大変複雑である。
noshi's 簡易 LARSCH 法は、LARSCH 法から "Reduce" Step を抜いたものである。

- LARSCH 法
  - メリット
    - 計算量が $O(N)$ とよい
  - デメリット
    - 実装が複雑
    - 辺コスト `f(i, j)` がスライドアクセスによってのみ高速に計算できる場合に対応できない
- noshi's 簡易 LARSCH 法
  - メリット
    - 実装が簡単
    - 辺コスト `f(i, j)` がスライドアクセスによってのみ高速に計算できる場合にも対応できる
  - デメリット
    - 計算量が $O(N \log N)$ と、LARSCH 法と比較すると若干悪い


## 仕様
構造体 `SimpleLARSCH<VAL>` を宣言して、関数 `solve(N, f)` を適用する。

- 関数 `solve()` の返り値: 受け取る変数を `res` として
  - `res[v].first`: 頂点 $0$ から頂点 $v$ への最短路の重み
  - `res[v].second`: 頂点 $0$ から頂点 $v$ への最短路において、$v$ に至る直前の頂点の番号 (存在しない場合は -1)
- また、構造体のメンバ変数にも結果が格納される
  - `dp[v]`: 頂点 $0$ から頂点 $v$ への最短路の重み
  - `cnt[v]`: 頂点 $0$ から頂点 $v$ への最短路に含まれる辺の本数
  - `prev[v]`: 頂点 $0$ から頂点 $v$ への最短路において、$v$ に至る直前の頂点の番号 (存在しない場合は -1)
- メンバ関数 `reconstruct()` を呼ぶことで、頂点 $0$ から頂点 $N$ への最短路が復元される


### ランダムアクセス版

```cpp
template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f, bool minimize_cnt = true)
```

- 整数 `N`: 頂点番号が $0, 1, 2, \dots, N$ であるとする (頂点数は $N + 1$)
- 関数 `f(l, r)`: 辺 $(l, r)$ のコストを返す関数


### スライドアクセス版

```cpp
template<class STATE, class ADD, class DEL, class GETCOST> vector<pair<VAL, int>> solve
(int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, bool minimize_cnt = true)
```

- テンプレート型 `STATE`: 辺 $(l, r)$ に対応する区間 $[l, r)$ についての情報をまとめた構造体の型
- 整数 `N`: 頂点番号が $0, 1, 2, \dots, N$ であるとする (頂点数は $N + 1$)
- `STATE` 型変数 `ini`: 幅のない区間に対応する値
- 関数 `add(st, i)`: 辺 $(l, r)$ のコストが、区間 $[l, r)$ についての値に対応するとして、区間に index $i$ の分を追加するときの `STATE` 型変数 `st` の変化を記述する関数
- 関数 `del(st, i)`: 辺 $(l, r)$ のコストが、区間 $[l, r)$ についての値に対応するとして、区間から index $i$ の分を削除するときの `STATE` 型変数 `st` の変化を記述する関数
- 関数 `get(st)`: 辺 $(l, r)$ のコストが、区間 $[l, r)$ についての値に対応するとして、区間の情報を表す `STATE` 型変数 `st` から辺コスト $f(l, r)$ を得るための関数


## 問題例

- [AtCoder EDPC Z - Frog 3](https://atcoder.jp/contests/dp/tasks/dp_z)
- [Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry](https://codeforces.com/contest/319/problem/C)
- [yukicoder No.705 ゴミ拾い Hard](https://yukicoder.me/problems/no/705)

## 参考文献

- tatyam: [Monge の手引き書](https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu)
- noshi: [簡易版 LARSCH Algorithm](https://noshi91.hatenablog.com/entry/2023/02/18/005856)
