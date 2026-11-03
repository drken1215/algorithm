---
title: UnionFind
documentation_of: ./union_find.hpp
---


## 概要

Union-Find
union by size による実装

## 使い方

- same(x, y): x と y が同じグループにいるか, 計算量はならし O(α(N))
- merge(x, y): x と y を同じグループにする, 計算量はならし O(α(N))
- size(x): x を含むグループの所属メンバー数, 計算量はならし O(α(N))
- group(x): x を含むグループを返す, 計算量は O(S) (S はグループのサイズ)
- groups(): グループ分けの構造を返す, 計算量は O(N)

