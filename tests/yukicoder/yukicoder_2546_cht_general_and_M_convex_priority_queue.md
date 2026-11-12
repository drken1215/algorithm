---
title: yukicoder No.2546 Many Arithmetic Sequences
documentation_of: ./yukicoder_2546_cht_general_and_M_convex_priority_queue.test.cpp
---

## 問題概要
$N$ 個の等差数列がある。等差数列 $i$ は初項 $A_i$、交差 $D_i$ である。

各等差数列 $i$ から先頭 $X_{i}$ 個を選ぶ。ここで、$X_{1} + X_{2} + \dots + X_{N} = M$ となるようにする。

選んだ数の総和の最大値を求めよ。


### 制約

- $N, M \le 3 \times 10^{5}$


## 解法メモ
めちゃくちゃ面白い！！

+ 交差が $0$ 以上であるものについては、1 箇所だけから選ぶようにするのが良いことが言える（「交換しても悪化しない」の議論）
+ 交差が負であるものについては、いわゆる資源配分問題（M 凹関数の和の最大化）になっており、priority_queue を用いた貪欲法が最適であることが言える

交差が $0$ 以上のところから合計 $0, 1, 2, \dots$ 個を選ぶ方法と、交差が負のところから合計 $0, 1, 2, \dots$ 個を選ぶ方法をそれぞれ求めてあげて、その結果をマージすればよい。
交差が負のところは簡単。交差が $0$ 以上のところを考える。


### 交差が $0$ 以上のところについて
関数

$f_{i}(x) = \frac{x}{2}(D_{i} x + 2A_{i} - D_{i})$

を考えて、$x = 0, 1, 2, \dots, m$ について $\min_{i}f_{i}(x)$ を求める問題である。いわゆる行最小値問題である。
$x$ を固定したときに、

$\min_{i} (D_{i} x + 2A_{i} - D_{i})$

を求めることができればよいことがわかる。これは一次関数群なので Convex Hull Trick (Li Chao Tree) で求めることができる。


### 関数群 $f_1, f_2, \dots$ の性質

- 関数群 $f_1, f_2, \dots$ は Monge ではない
- 関数群 $f_1, f_2, \dots$ は定義域全体で見ると Pairwise Single Crossing でもない
- しかし、$x > 0$ の領域に限定すると、関数群 $f_1, f_2, \dots$ は Pairwise Single Crossing である
  - なので、一次関数群に帰着しなくても General Li Chao Tree で解ける


