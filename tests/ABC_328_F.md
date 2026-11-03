---
title: 列挙可能 Union-Find の例 (ABC 328 F)
documentation_of: ./ABC_328_F.cpp
---

## 概要

Union-Find

- union by size + 経路圧縮による実装
- `nex` を用いて、x を含むグループを取り出す関数 `group(x)` を備える (通称：列挙可能 Union-Find)

## 仕様

- `same(x, y)`: x と y が同じグループにいるか, 計算量はならし O(α(N))
- `merge(x, y)`: x と y を同じグループにする, 計算量はならし O(α(N))
- `size(x)`: x を含むグループの所属メンバー数, 計算量はならし O(α(N))
- `group(x)`: x を含むグループを返す, 計算量は O(S) (S はグループのサイズ)
- `groups()`: グループ分けの構造を返す, 計算量は O(N)

## 参考文献

- noshi: [素集合の要素の列挙と併合 (単方向循環リスト)](https://noshi91.hatenablog.com/entry/2019/07/19/180606)
