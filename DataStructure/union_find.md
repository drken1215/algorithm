---
title: union_find
documentation_of: ./union_find.hpp
---

## 概要

Union-Find

- union by size + 経路圧縮による実装
- `nex` を用いて、$x$ を含むグループを取り出す関数 `group(x)` を備える (通称：列挙可能 Union-Find)

## 仕様

- `same(x, y)`: $x$ と $y$ が同じグループにいるか。計算量はならし $O(\alpha(N))$
- `merge(x, y)`: $x$ と $y$ を同じグループにする。計算量はならし $O(\alpha(N))$
- `size(x)`: $x$ を含むグループの所属メンバー数。計算量はならし $O(\alpha(N))$
- `group(x)`: $x$ を含むグループを返す。計算量は $O(S)$ ($S$ はグループのサイズ)
- `groups()`: グループ全体の構造を返す。計算量は $O(N \alpha(N))$

## 問題例

- [Library Checker - Unionfind](https://judge.yosupo.jp/problem/unionfind)
- [AtCoder Library Practice Contest A - Disjoint Set Union](https://atcoder.jp/contests/practice2/tasks/practice2_a)
- [AtCoder PAST 第10回 H - 連結成分](https://atcoder.jp/contests/past202203-open/tasks/past202203_h) (`group(x)` を使用)
- [AtCoder ABC 328 F - Good Set Query](https://atcoder.jp/contests/abc328/tasks/abc328_f) (`group(x)` を使用)

## 参考文献

- noshi: [素集合の要素の列挙と併合 (単方向循環リスト)](https://noshi91.hatenablog.com/entry/2019/07/19/180606)
