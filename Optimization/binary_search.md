---
title: binary_search
documentation_of: ./binary_search.hpp
---

## 概要

めぐる式二分探索。

|ok - ng| = 1 となるような (ng, ok) を見つけて、そのときの ok を返す

- `binary_search_01(ng, ok, f)`
- `binary_search_10(ok, ng, f)`

## 制約

- `f(ng)` = False
- `f(ok)` = True

## 問題例

- [競プロ典型 90 問 001 - Yokan Party（★4）](https://atcoder.jp/contests/typical90/tasks/typical90_a)

## 参考文献

- drken: [二分探索アルゴリズムを一般化 〜 めぐる式二分探索法のススメ 〜](https://qiita.com/drken/items/97e37dd6143e33a64c8c)
