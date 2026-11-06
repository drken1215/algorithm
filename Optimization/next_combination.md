---
title: next_combination (nCk 通りの全探索)
documentation_of: ./next_combination.hpp
---

## 概要
台集合 $V = ${$1, 2, \dots, N$} に対して、サイズの等しい $V$ の部分集合を列挙したいとする。

$V$ の部分集合 $S$ に対して、$S$ とサイズの等しい部分集合を辞書順に列挙したときの、「$S$ の次の部分集合」を求める。

## 仕様

- `next_combination(S, N)`: サイズ $N$ の台集合 $V$ の部分集合 $S$ を、$S$ とサイズの等しい辞書順で次のものに置き換える。
  - 次のものが存在するならば `True` を返し、存在しない場合は `False` を返す。

## 使い方
たとえば $N = 5$ として、サイズ $3$ の部分集合を列挙したいとき、次のように書ける。

```cpp
#include <bits/stdc++.h>
using namespace std;

#include "bit_utility.hpp"

int main() {
    int N = 5;
    int S = 0b00000111;
    do {
        cout << bitset<5>(S) << endl;
    } while (next_combination(S, N));
}
```

### 実行結果

```
00111
01011
01101
01110
10011
10101
10110
11001
11010
11100
```

## 問題例

- [AtCoder ABC 328 E - Modulo MST](https://atcoder.jp/contests/abc328/tasks/abc328_e)

## 参考文献

