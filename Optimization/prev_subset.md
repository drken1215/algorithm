---
title: prev_subset
documentation_of: ./prev_subset.hpp
---

## 概要
台集合 $V = \{1, 2, \dots, N\}$ の部分集合 $U$ に対して、$U$ の部分集合を列挙したいとする。

集合 $U$ の部分集合を辞書順逆順に列挙したときの、「部分集合 $S$ の前の部分集合」を求める。

## 仕様

- `prev_subset(S, U)`: 集合 $U$ の部分集合 $S$ を、辞書順で前のものに置き換える。
  - 前のものが存在するならば `True` を返し、存在しない場合は `False` を返す。

## 問題例

## 参考文献

- drken: [ビット演算 (bit 演算) の使い方を総特集！ 〜 マスクビットから bit DP まで 〜](https://qiita.com/drken/items/7c6ff2aa4d8fce1c9361)
