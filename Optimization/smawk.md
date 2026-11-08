---
title: TM 行最小値問題 (SMAWK 法, in O(H + W))
documentation_of: ./smawk.hpp
---

## 概要
$H \times W$ 行列 $A$ が Totally Monotone であるとき、
計算量 $O(H + W)$ で、$A$ の各行に対する最小値を求める。

## 仕様

```cpp
template<class VAL, class FUNC> vector<pair<VAL, int>> SMAWK(int H, int W, const FUNC &f)
```

- $H \times W$ 行列の $(i, j)$ 成分が $f(i, j)$ で表されるとする
- 返り値: 受け取る変数を `res` として
  - `res[l].first`: $l$ 行目の最小値
  - `res[l].second`: $l$ 行目の値が最小となる列の番号


## 問題例

- [COLOCON 2018 Final C - スペースエクスプローラー高橋君](https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c)

## 参考文献

- tatyam: [Monge の手引き書](https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu)
