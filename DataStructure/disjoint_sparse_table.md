---
title: Disjoint Sparse Table
documentation_of: ./disjoint_sparse_table.hpp
---

## 概要
Sparse Table は本来は「交叉半束」上で動作するが、Disjoint Sparse Table は「半群」上でも動作できるようにしたものである。
 
動的更新はできないが、最初に静的に構築してしまえば、区間取得を $O(1)$ で実行できる。

## 仕様

- `init(vec)`: 初期化構築, $O(N \log N)$
- `get(a, b)`: 区間 $[a, b)$ の最小値を取得する, $O(1)$

## 制約

- 半群: 二項演算が定義されて、結合法則を満たす、単位元も逆元も不要 (min, max, 和, 積, gcd)
  - モノイド: 半群 + 単位元
  - 群: 半群 + 単位元 + 逆元
- 交叉半束: さらに「べき等則 (a + a = a)」「可換則 (a + b = b + a)」を満たす必要がある

## 問題例

- [Codechef Product on the segment by modulo](https://www.codechef.com/problems/SEGPROD)

## 参考文献

- noshi: [Disjoint Sparse Table と セグ木に関するポエム](http://noshi91.hatenablog.com/entry/2018/05/08/183946)
