---
title: yukicoder No.2404 Vertical Throw Up
documentation_of: ./yukicoder_2404_general_li_chao_tree.test.cpp
---

## 問題概要
$a$ は正の定数とする。以下の $Q$ 個のクエリに答えよ。

- クエリタイプ 1 ($s, t$)：放物線 $y = -a(x - s)(x - t)$ が新たに挿入される
- クエリタイプ 2 ($t$)：すでに挿入された各放物線について、$x = t$ における値を求め、その最小値を答えよ

### 制約

- $2 \le Q \le 10^{5}$
- $0 \le s, t \le 10^{8}$


## 解法メモ
クエリタイプ 1 で挿入される関数を、挿入順に番号をつけて $f_1(x), f_2(x), \dots$ とする。

- 関数群 $f_1(x), f_2(x), \dots$ は Monge ではない
- 関数群 $f_1(x), f_2(x), \dots$ は Pairwise Single Crossing ではある

よって、Convex Hull Trick (General Li Chao Tree) で解ける。

