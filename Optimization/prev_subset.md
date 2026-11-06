---
title: prev_subset
documentation_of: ./prev_subset.hpp
---

## 概要
台集合 $V = ${$1, 2, \dots, N$} の部分集合 $U$ に対して、$U$ の部分集合を列挙したいとする。

集合 $U$ の部分集合を辞書順逆順に列挙したときの、「部分集合 $S$ の前の部分集合」を求める。

## 仕様

- `prev_subset(S, U)`: 集合 $U$ の部分集合 $S$ を、辞書順で前のものに置き換える。
  - 前のものが存在するならば `True` を返し、存在しない場合は `False` を返す。

## 使い方
たとえば $N = 5$ として、集合 $S$ = {$0, 2, 4$} の部分集合を列挙したいとき、次のように書ける。

```cpp
#include <bits/stdc++.h>
using namespace std;

#include "bit_utility.hpp"

int main() {
    int S = 0b10101;
    int v = S;
    do {
        cout << bitset<5>(v) << endl;
    } while (prev_subset(v, S));
}
```

### 実行結果

```
10101
10100
10001
10000
00101
00100
00001
00000
```

## 問題例

## 参考文献

- drken: [ビット演算 (bit 演算) の使い方を総特集！ 〜 マスクビットから bit DP まで 〜](https://qiita.com/drken/items/7c6ff2aa4d8fce1c9361)
