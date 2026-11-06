---
title: 数独ソルバー
documentation_of: ./sudoku.cpp
---

## 概要
数独の問題を標準入力として与えるとそれを解いて出力します。

数独は、3 × 3 ブロックに区切られた 9 × 9 の正方形の各マスに 1 から 9 までの数値を入れるペンシルパズルです。縦・横の各列および、太線で囲まれた 3 × 3 のブロック内には、同じ数値が複数入らないようにします。

<img src=https://github.com/drken1215/sudoku/blob/master/sudoku_example.png width=500mm>


## 仕様

### 入力
以下のいずれかをコンストラクタで与える

- vector<string> 型変数 (空きマスはデフォルトでは '*')
- vector<vector<int>> 型変数 (空きマスはデフォルトでは '-1')

### 使用例

```cpp
vector<string> input = {
    "53**7****",
    "6**195***",
    "*98****6*",
    "8***6***3",
    "4**8*3**1",
    "7***2***6",
    "*6****28*",
    "***419**5",
    "****8**79"
};
Sudoku sudoku(input);

// 数独を解く (すべての解を求める)
vector<vector<vector<int>>> res = sudoku.solve();
```

## 問題例

- [POJ 2676 - Sudoku (solve 9 x 9 sudoku)](http://poj.org/problem?id=2676)
- [AtCoder ABC 327 C - Number Place (judge validity)](https://atcoder.jp/contests/abc327/tasks/abc327_c)
