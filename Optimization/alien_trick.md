---
title: Alien Trick --- min f(x) s.t. g(x) = K (by ラグランジュ緩和)
documentation_of: ./alien_trick.hpp
---

## 概要

### 問題設定
変数 $x$ について、なんらかの「個数」を表す関数を $g(x)$ とおく（たとえば $x$ がグラフのパスを表し、$g(x)$ がパス中の辺の本数を表すなど）。
$g(x)$ は $L$ 以上 $R$ 以下の整数値をとる関数とする。

$K$ を $L \le K \le R$ を満たす正の整数として、次の最適化問題を考える。

$\mathrm{min}_{x} f(x)$ $\mathrm{s.t.}$ $g(x) = K$

$g(x) = p$ となる $x$ についての $f(x)$ の最小値を $h(p)$ とする。求めたい値は $h(K)$ である。ここで、次の仮定をおく。

----
$h(p)$ が $p$ について下に凸である

----

このとき、ラグランジュ強双対性が成り立つことが知られている（下の参考文献を参照）。
すなわち、ラグランジュ双対変数 $\lambda$ を整数として、次のラグランジュ緩和問題を考えると、この最適値はもとの最適化問題の最適値に一致する。

$\mathrm{max}_{\lambda}$ $\mathrm{min}_{x}$ $(f(x) + \lambda(g(x) - K))$


### ラグランジュ緩和問題の解き方
$\lambda$ を固定して、最適化問題

$P(\lambda): \mathrm{min}_{x} (f(x) + \lambda g(x))$

を考える（最適値を $L(\lambda)$ とおく）。
この最適化問題 $P(\lambda)$ の最適解の集合は、$\mathrm{min}_{x} (f(x) + \lambda(g(x) - K))$ の最適解の集合に一致する。

最適化問題 $P(\lambda)$ の最適解 $x^{\mathrm{opt}}$ を 1 つとり、$k(\lambda) = g(x^{\mathrm{opt}})$ とおくと、
$k(\lambda)$ は $\lambda$ について広義単調減少となることが知られている（[Alien Trick](https://info.atcoder.jp/entry/algorithm_lectures/alien_trick) を参照）。

よって、二分探索法によって、

$k(\lambda_{1}) \le K \le k(\lambda_{2})$ ($\lambda_{1} = \lambda_{2} + 1$)

を満たすような $\lambda^{*}$ を求めることができる。
このとき、$\lambda = \lambda_{1}, \lambda_{2}$ のいずれかについて、
$g(x) = K$ であるような最適化問題 $P(\lambda)$ の最適解 $x$ が存在することが示せる。
すなわち、

$h(K) = \min(L(\lambda_{1}) - \lambda_{1} K, L(\lambda_{2}) - \lambda_{2} K)$

が成立する。


## 仕様

```cpp
template<class FUNC> pair<VAL, VAL> solve(const FUNC &lag, VAL K)
```

- `lag(lambda)`: $x^{\mathrm{opt}} = \mathrm{argmin}_{x} (f(x) + \lambda g(x))$ としたときのペア値 $(f(x^{\mathrm{opt}}), g(x^{\mathrm{opt}}))$ を返す関数
- 返り値: (元の最適化問題の最適値, そのときの $\lambda$)


## 問題例

- [ABC 218 H - Red and Blue Lamps](https://atcoder.jp/contests/abc218/tasks/abc218_h)
- [Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem](https://codeforces.com/contest/868/problem/F)

## 参考文献

- AtCoder Algorithm Lectures: [Alien Trick](https://info.atcoder.jp/entry/algorithm_lectures/alien_trick)
- noshi: [Monge グラフ上の d-辺最短路長を計算するアルゴリズム](https://noshi91.github.io/algorithm-encyclopedia/d-edge-shortest-path-monge)
- noshi: [Aliens DP における二分探索の色々](https://noshi91.hatenablog.com/entry/2023/11/20/052227)
- noshi: [Aliens DP で辺の本数が区間で指定される場合](https://noshi91.hatenablog.com/entry/2022/01/13/001217)
- wata: [双対性](https://www.slideshare.net/slideshow/ss-91375739/91375739#1)
