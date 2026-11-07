---
title: Sparse Table (スパーステーブル)
documentation_of: ./sparse_table.hpp
---

## 概要
動的更新はできないが、最初に静的に構築してしまえば、区間最小値を $O(1)$ でとり出すことができる。

## 仕様

- `init(vec)`: 初期化構築, $O(N \log N)$
- `get(a, b)`: 区間 $[a, b)$ の最小値を取得する, $O(1)$

## 制約

- 交叉半束 (meet-semi-lattice) 上で動作する
- 半群 (semi-group) 上で動作する上位互換バージョン (disjoint sparse table) もある

## 問題例

- [AtCoder ABC 407 F - Sums of Sliding Window Maximum](https://atcoder.jp/contests/abc407/tasks/abc407_f)

## 参考文献

