---
title: General Li Chao Tree (一般関数対応)
documentation_of: ./general_li_chao_tree.hpp
---

## 概要
一般関数群 (Pairwise Single Crossing) に対応した Li Chao Tree。

- Dynamic 版 ($x$ 座標を座標圧縮せず、有効な定義域を与える。$x$ 座標は $10^{9}$ や $-10^{9}$ などでも構わない)
- 座標圧縮版 (クエリを呼び出す $x$ 座標を座標圧縮する)

を用意している。後者はクエリ先読みができることを前提としているが、基本的には前者よりも高速である。


## 制約

- 挿入する関数群 $f_1, f_2, \dots, f_N$ が Pairwise Single Crossing であること
  - 「$f_1, f_2, \dots, f_N$ を適切に並び替えると **Monge** になること」は十分条件である
  - 「$f_1, f_2, \dots, f_N$ を適切に並び替えると **anti-Monge** になること」も十分条件である


## 仕様 (Dynamic 版)

```cpp
template<class VAL, class FUNC, class X = long long> struct DynamicGeneralLiChaoTree;
```

- テンプレート引数:
  - `VAL`: 関数の $f(x)$ の値の型
  - `FUNC`: 関数の型 (具体的な指定方法は使用例を参照)
  - `X`: $x$ 座標の型
- コンストラクタ: `DynamicGeneralLiChaoTree(X xl, X xr, VAL inf = numeric_limits<VAL>::max() / 2)`
  - `[xl, xr)`: $x$ の範囲 (この幅を $X$ とする)
  - `inf`: 無限大を表す値
- 関数の挿入
  - `insert(f)`: 関数 $f$ を挿入する。計算量 $O(\log X)$ 
  - `insert_segment(f, l, r)`: 関数 $f$ を区間 $[l, r)$ の範囲内で挿入する。計算量 $O(\log^{2} X)$
- 1 点取得
  - `query(x)`: {$\min_{i} f_{i}(x)$, 最小を達成する $i$} を返す ($x$ に対して定義された関数が存在しない場合は {INF, -1})。計算量 $O(\log X)$


### Dynamic 版の使用例
以下のように、関数群を「ラムダ式関数を返すファクトリ」として定義する。

```cpp
// 関数群の定義
auto make_func = [&](long long a, long long b) {
    return [&, a, b](long long x) -> long long {
        return a * x + b;
    };
};

// Li Chao Tree の定義
long long LX = -(1LL << 30), RX = 1LL << 30;
DynamicGeneralLiChaoTree<long long, decltype(make_func(0, 0))> cht(LX, RX);

// 区間挿入
cht.insert_segment(make_func(a, b), l, r);

// クエリ取得
auto [res, id] = cht.query(p);
```


## 仕様 (座標圧縮版)

```cpp
template<class VAL, class FUNC, class X = long long> struct GeneralLiChaoTree;
```

- テンプレート引数:
  - `VAL`: 関数の $f(x)$ の値の型
  - `FUNC`: 関数の型 (具体的な指定方法は使用例を参照)
  - `X`: $x$ 座標の型
- コンストラクタ 1: `GeneralLiChaoTree(const vector<X>& xs, VAL inf = numeric_limits<VAL>::max() / 2, bool sorted = true)`
  - `xs`: `query(x)` の対象となる $x$ を重複なく小さい順に格納したもの (この値の個数を $N$ とする)
    - 区間挿入 `insert_segment(f, l, r)` の端点 $l, r$ については、`xs` に含まれている必要はない
  - `inf`: 無限大を表す値
  - `sorted`: `xs` がソート済みであるかどうかを表すフラグ変数
- コンストラクタ 2: `GeneralLiChaoTree(int xl, int xr, VAL inf = numeric_limits<VAL>::max() / 2)`
  - `[xl, xr)`: $x$ の範囲 ($N = \mathrm{xr} - \mathrm{xl}$ とする) 
  - `inf`: 無限大を表す値
- 関数の挿入
  - `insert(f)`: 関数 $f$ を挿入する。計算量 $O(\log N)$ 
  - `insert_segment(f, l, r)`: 関数 $f$ を区間 $[l, r)$ の範囲内で挿入する。計算量 $O(\log^{2} N)$
- 1 点取得
  - `query(x)`: {$\min_{i} f_{i}(x)$, 最小を達成する $i$} を返す ($x$ に対して定義された関数が存在しない場合は {INF, -1})。計算量 $O(\log N)$


### 座標圧縮版の使用例

```cpp
// 関数群の定義
auto make_func = [&](long long a, long long b) {
    return [&, a, b](long long x) -> long long {
        return a * x + b;
    };
};

// Li Chao Tree の定義
vector<int> xs;  // query(x) を呼ぶ x の集合 (sort されていて、unique になっている必要がある)
DynamicGeneralLiChaoTree<long long, decltype(make_func(0, 0))> cht(xs);

// 区間挿入
cht.insert_segment(make_func(a, b), l, r);

// クエリ取得
auto [res, id] = cht.query(p);
```


## 問題例

- [Library Checker - Line Add Get Min](https://judge.yosupo.jp/problem/line_add_get_min) (直線を管理)
- [Library Checker - Segment Add Get Min](https://judge.yosupo.jp/problem/segment_add_get_min) (直線を管理。区間挿入あり)
- [COLOCON 2018 Final C - スペースエクスプローラー高橋君](https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c) (一般関数)
- [AtCoder EDPC Z - Frog 3](https://atcoder.jp/contests/dp/tasks/dp_z) (Monge 単一始点最短路)
- [Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry](https://codeforces.com/contest/319/problem/C) (Monge 単一始点最短路)
- [yukicoder No.705 ゴミ拾い Hard](https://yukicoder.me/problems/no/705) (Monge 単一始点最短路)
- [AtCoder ABC 480 G - Zonal Score Minimization](https://atcoder.jp/contests/abc480/tasks/abc480_g) (anti-Monge 単一始点最短路。定義域全体に anti-Monge 性を拡張する方針でも、区間挿入を使っても解ける)


## 参考文献

- Algorithms for Competitive Programming: [Convex hull trick and Li Chao tree](https://cp-algorithms.com/geometry/convex_hull_trick.html)
- robert1003: [A Simple Introduction to Li-Chao Segment Tree](https://robert1003.github.io/2020/02/06/li-chao-segment-tree.html)
- rama_pang: [[Tutorial] Li Chao Tree Extended](https://codeforces.com/blog/entry/86731)
