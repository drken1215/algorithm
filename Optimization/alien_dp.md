---
title: Monge 単一始点 d-辺最短路問題 (by Alien DP)
documentation_of: ./alien_dp.hpp
---

## 概要
頂点番号が $0, 1, 2, \dots, N$ であるような頂点数 $N$ の DAG を考える。
$0 \le l < r \le N$ を満たす $l, r$ について、有向辺 $(l, r)$ があり、重み $f(l, r)$ は Monge であるとする。

このとき、頂点 $0$ から頂点 $N$ へと至る、辺数 $D$ 本のパスの重みの最小値を [Alien Trick](https://drken1215.github.io/algorithm/Optimization/alien_trick.hpp.html) を用いて解く。
この解法をとくに Alien DP と呼ぶ。

Monge 単一始点最短路問題の計算量を $O(f(N))$ とすると、Alien DP の計算量は $O(f(N) \log C)$ である ($C$ は辺の重みの絶対値の最大値)。


## 仕様

### ランダムアクセス版

```cpp
template<class FUNC> pair<VAL, VAL> solve(int N, const FUNC &f, int D, const string &solver = "simple_larsch")
```

- 整数 `N`: 頂点番号が $0, 1, 2, \dots, N$ であるとする (頂点数は $N + 1$)
- 関数 `f(l, r)`: 辺 $(l, r)$ のコストを返す関数
- 整数 `D`: 含まれる辺数がちょうど $D$ 本であるパスを考える
- 文字列 `solver`: "monotone", "smawk", "larsch", "simple_larsch" のいずれかを指定する
- 返り値: (最短路長, ラグランジュ変数 $\lambda$ の値)


### スライドアクセス版

```cpp
template<class STATE, class ADD, class DEL, class GETCOST> pair<VAL, int> solve
(int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, int D, const string &solver = "simple_larsch")
```

- 整数 `N`: 頂点番号が $0, 1, 2, \dots, N$ であるとする (頂点数は $N + 1$)
- `STATE` 型変数 `ini`: 幅のない区間に対応する値
- 関数 `add(st, i)`: 辺 $(l, r)$ のコストが、区間 $[l, r)$ についての値に対応するとして、区間に index $i$ の分を追加するときの `STATE` 型変数 `st` の変化を記述する関数
- 関数 `del(st, i)`: 辺 $(l, r)$ のコストが、区間 $[l, r)$ についての値に対応するとして、区間から index $i$ の分を削除するときの `STATE` 型変数 `st` の変化を記述する関数
- 関数 `get(st)`: 辺 $(l, r)$ のコストが、区間 $[l, r)$ についての値に対応するとして、区間の情報を表す `STATE` 型変数 `st` から辺コスト $f(l, r)$ を得るための関数
- 整数 `D`: 含まれる辺数がちょうど $D$ 本であるパスを考える
- 文字列 `solver`: "monotone", "smawk", "larsch", "simple_larsch" のいずれかを指定する
- 返り値: (最短路長, ラグランジュ変数 $\lambda$ の値)


## 問題例

- [ABC 218 H - Red and Blue Lamps](https://atcoder.jp/contests/abc218/tasks/abc218_h)
- [Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem](https://codeforces.com/contest/868/problem/F)

## 参考文献

- AtCoder Algorithm Lectures: [Alien Trick](https://info.atcoder.jp/entry/algorithm_lectures/alien_trick)
- noshi: [Monge グラフ上の d-辺最短路長を計算するアルゴリズム](https://noshi91.github.io/algorithm-encyclopedia/d-edge-shortest-path-monge)
- noshi: [Aliens DP における二分探索の色々](https://noshi91.hatenablog.com/entry/2023/11/20/052227)
- noshi: [Aliens DP で辺の本数が区間で指定される場合](https://noshi91.hatenablog.com/entry/2022/01/13/001217)
- wata: [双対性](https://www.slideshare.net/slideshow/ss-91375739/91375739#1)
