---
title: 列挙可能 Union-Find の例 (ABC 328 F)
documentation_of: ./ABC_328_F_enumerable_union_find.cpp
---

## 列挙可能 Union-Find の概要

「列挙可能 Union-Find」の使用例である。

"[union_find.hpp](https://drken1215.github.io/algorithm/DataStructure/union_find.hpp.html)" に列挙機能を搭載している。

- union by size + 経路圧縮による実装
- `nex` を用いて、`x` を含むグループを取り出す関数 `group(x)` を備える (通称：列挙可能 Union-Find)


## ABC 328 F の概要
$𝑁$ 個の整数値 $𝑥_1,𝑥/2, \dots, 𝑥_𝑁$ に関する制約条件が $𝑄$ 個与えられる。

$𝑖$ 番目の制約条件では 3 つの整数の組 $(𝑎,𝑏,𝑑)$ が与えられ、$𝑥_𝑎 − 𝑥_𝑏 = 𝑑$ という形をしている。ここで、次のクエリに答えよ。

すでにある制約条件と、この $𝑖$ 番目の制約条件「$𝑥_𝑎 − 𝑥_𝑏 = 𝑑$」とが矛盾するかどうかを判定せよ。
矛盾しない場合には $𝑖$ を出力し、矛盾する場合には何も出力しないこととする。


## ABC 328 F の解法
制約条件 $𝑥_𝑎 − 𝑥_𝑏 = 𝑑$ によって要素 $𝑎,𝑏$ をマージする際に、根付き木のサイズが小さい方の各頂点について、$𝑥$ の値を調整する。

マージテクの仕組みによって、全体の計算量は $𝑂⁢(𝑁 + 𝑄⁢ \log⁡ 𝑁)$ となる。

### 別解
ポテンシャル付き Union-Find を用いる解法もある


## 参考文献

- けんちょん: [AtCoder ABC 328 F - Good Set Query (1D, 水色, 525 点)](https://drken1215.hatenablog.com/entry/2023/11/12/040700)
- noshi: [素集合の要素の列挙と併合 (単方向循環リスト)](https://noshi91.hatenablog.com/entry/2019/07/19/180606)
