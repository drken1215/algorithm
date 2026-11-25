---
title: dynamic_modint
documentation_of: ./dynamic_modint.hpp
---

## 概要
通常の modint は、法を表す素数 $p$ が定数であることを想定している ($998244353$ や $1000000007$ など)。これに対して、法を表す素数 $p$ が入力で与えられる場合がある。

このような場合に対応した、有限体 $\mathbb{F}_p$ を表す構造体である。通常の modint よりも汎用性は高いが、実行速度がやや劣る。
