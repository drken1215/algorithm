---
title: Alien Trick --- min f(x) s.t. g(x) = K (by ラグランジュ緩和)
documentation_of: ./alien_trick.hpp
---

## 概要

### 問題設定
変数 $x$ について、なんらかの「個数」を表す関数を $g(x)$ とおく（たとえば $x$ がグラフのパスを表し、$g(x)$ がパス中の辺の本数を表すなど）。
$g(x)$ は $L$ 以上 $R$ 以下の整数値をすべてとる関数とする。

$K$ を $L \le K \le R$ を満たす正の整数として、次の最適化問題を考える ($f(x)$ は整数値をとる関数とする)。

$\min_{x} f(x)$ $\mathrm{s.t.}$ $g(x) = K$

$g(x) = k$ となる $x$ についての $f(x)$ の最小値を $h(k)$ とする。
$h(k)$ の定義域は $k = L, L+1, \dots, R$ であり、求めたい値は $h(K)$ である。
ここで、次の仮定をおく。

----

$h(k)$ が定義域内で下に凸であり、離散凸性を満たす。すなわち、

$h(k-1) + h(k+1) \ge 2h(k)$

を満たす。

----

この仮定を満たす具体的な状況については、[Alien Trick](https://info.atcoder.jp/entry/algorithm_lectures/alien_trick) で、

- ソートして交互にとる
- ネットワークフローとマッチング
- Monge 性 (とくに、Alien DP と呼ばれる解法に繋がる)
- マトロイド

という 4 つの場面が紹介されている。


### ラグランジュ緩和した問題
上記の仮定の下で、ラグランジュ強双対性が成り立つことが知られている。

- [Monge グラフ上の d-辺最短路長を計算するアルゴリズム](https://noshi91.github.io/algorithm-encyclopedia/d-edge-shortest-path-monge) では、ラグランジュ双対性の言葉を意識した証明
- [Alien Trick](https://info.atcoder.jp/entry/algorithm_lectures/alien_trick) では、凸数列の凸共役性の言葉を意識した証明

を与えている。
具体的には、ラグランジュ双対変数 $\lambda$ を整数として、ラグランジュ緩和した問題

$\max_{\lambda} \min_{x} (f(x) + \lambda(g(x) - K))$

の最適値は、もとの最適化問題の最適値に一致する。補足として、ラグランジュ強双対性が成り立たなくても、以下のことが成り立つ。

- 一般に、ラグランジュ緩和した問題の最適解の値は、もとの最適化問題の下界を与える (弱双対性)
- 一般に、$\min_{x} (f(x) + \lambda(g(x) - K))$ は $\lambda$ について上に凸であり、ラグランジュ緩和した問題の最適解は三分探索によって求められる

ここでは、三分探索法よりも高速な二分探索法を考える。


### ラグランジュ緩和した問題の解き方
$\lambda$ を固定して、最適化問題

$P(\lambda): \min_{x} (f(x) + \lambda g(x))$

を考える（最適値を $L(\lambda)$ とおく）。
この最適化問題 $P(\lambda)$ は、$\min_{x} (f(x) + \lambda(g(x) - K))$ と等価であることに注意する。
また、最適化問題 $P(\lambda)$ には、Fenchel 双共役性に基づく幾何学的解釈も存在する（[Alien Trick](https://info.atcoder.jp/entry/algorithm_lectures/alien_trick) を参照）。

さて、$g(x) = K$ となるように $\lambda$ をコントロールすることを考える。
最適化問題 $P(\lambda)$ の最適解 $x^{\mathrm{opt}}$ を 1 つとり、$k(\lambda) = g(x^{\mathrm{opt}})$ とおくと、
$k(\lambda)$ は $\lambda$ について広義単調減少となる（[Alien Trick](https://info.atcoder.jp/entry/algorithm_lectures/alien_trick) を参照）。

よって、二分探索法によって、

$k(\lambda_{1}) \le K \le k(\lambda_{2})$ ($\lambda_{1} = \lambda_{2} + 1$)

を満たすような $\lambda_{1}, \lambda_{2}$ を求めることができる。
このとき、$\lambda = \lambda_{1}, \lambda_{2}$ のいずれかについて、
$g(x) = K$ であるような最適化問題 $P(\lambda)$ の最適解 $x$ が存在することが示せる。
すなわち、

$h(K) = \max(L(\lambda_{1}) - \lambda_{1} K, L(\lambda_{2}) - \lambda_{2} K)$

が成立する。


### 復元について
$k(\lambda_{1}) \le K \le k(\lambda_{2})$ ($\lambda_{1} = \lambda_{2} + 1$)

を満たすような $\lambda_{1}, \lambda_{2}$ を求めることができても、実際に

$\min_{x} f(x)$ $\mathrm{s.t.}$ $g(x) = K$

の最適解 $x$ を復元することは、一般には難しいと考えられる。ただし、以下に挙げるような、最適解に交換性があるような個別の場合については、最適解を復元できる場合がある。

- 最小費用流問題で、流量を 1 ずつ増減できる場合 (Alien Trick ではないが、[ABC 373 G - Unevenness](https://atcoder.jp/contests/abc393/tasks/abc393_g) など)
- マトロイドで、基の交換性を利用できる場合
- DP の最適遷移を組み合わせて、選択肢を 1 ずつ変えられる場合

これらの場合には、$g(x)$ が $K$ よりも大きい最適解と小さい最適解とから、最適性を保ちながら $g(x) = K$ であるような解を構成できる。


## 仕様

```cpp
template<class VAL> struct AlienTrick;
```

- `VAL` は符号付き整数型とする
- オーバーフローに注意すること


### 制約条件が g(x) = K の場合

```cpp
template<class FUNC> pair<VAL, VAL> solve(const FUNC &lag, VAL K)
```

- `lag(lambda)`: $x^{\mathrm{opt}} = \mathrm{argmin}_{x} (f(x) + \lambda g(x))$ としたときのペア値 $(f(x^{\mathrm{opt}}), g(x^{\mathrm{opt}}))$ を返す関数
  - 最大化問題の場合は $x^{\mathrm{opt}} = \mathrm{argmax}_{x} (f(x) + \lambda g(x))$ としたときのペア値 $(f(x^{\mathrm{opt}}), g(x^{\mathrm{opt}}))$ を返す関数
- 返り値: (元の最適化問題の最適値, そのときの $\lambda$)


### 制約条件が L <= g(x) <= R の場合

```cpp
template<class FUNC> pair<VAL, VAL> solve(const FUNC &lag, VAL L, VAL R) 
```

インターフェースについては同上。


## 使用例
Alien Trick は、とくに Monge 単一始点 d-辺最短路問題への応用 (Alien DP) がよく問われる。以下のページに、Alien DP における使用例を載せている。

- [Monge 単一始点 d-辺最短路問題 (by Alien DP)](https://drken1215.github.io/algorithm/Optimization/alien_dp.hpp.html)


## 問題例

- [ABC 218 H - Red and Blue Lamps](https://atcoder.jp/contests/abc218/tasks/abc218_h)
- [Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem](https://codeforces.com/contest/868/problem/F)

## 参考文献

- AtCoder Algorithm Lectures: [Alien Trick](https://info.atcoder.jp/entry/algorithm_lectures/alien_trick)
- noshi: [Monge グラフ上の d-辺最短路長を計算するアルゴリズム](https://noshi91.github.io/algorithm-encyclopedia/d-edge-shortest-path-monge)
- noshi: [Aliens DP における二分探索の色々](https://noshi91.hatenablog.com/entry/2023/11/20/052227)
- noshi: [Aliens DP で辺の本数が区間で指定される場合](https://noshi91.hatenablog.com/entry/2022/01/13/001217)
- wata: [双対性](https://www.slideshare.net/slideshow/ss-91375739/91375739#1)
