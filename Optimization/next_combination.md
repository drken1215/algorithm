---
title: next_combination
documentation_of: ./next_combination.hpp
---

## 概要
台集合 $V = \{1, 2, \dots, N\}$ に対して、サイズの等しい $V$ の部分集合を列挙したいとする。

$V$ の部分集合 $S$ に対して、$S$ とサイズの等しい部分集合を辞書順に列挙したときの、「$S$ の次の部分集合」を求める。

集合 $U$ の部分集合 $S$ を辞書順逆順に列挙したときの、「部分集合 $S$ の部分集合」を求める。

## 仕様

- `next_combination(S, N)`: サイズ $N$ の台集合 $V$ の部分集合 $S$ を、$S$ とサイズの等しい辞書順で次のものに置き換える。
  - 次のものが存在するならば `True` を返し、存在しない場合は `False` を返す。

## 問題例

- [AtCoder ABC 328 E - Modulo MST](https://atcoder.jp/contests/abc328/tasks/abc328_e)

## 参考文献

