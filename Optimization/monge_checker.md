---
title: Monge 性・anti-Monge 性の判定 (in O(HW))
documentation_of: ./monge_checker.hpp
---

## 概要
Monge 性・anti-Monge 性を判定する。

## 仕様

```cpp
template<class FUNC> bool is_monge
(int N, int M, const FUNC &f, bool upper_triangular = false, bool output_detail = true)
```

- 関数 $f(l, r)$ は、$0 \le l < N$, $0 \le r < M$ に対して定義されるとする
  - `upper_triangular` = `true` のときは、$0 \le l < r \le N$ を満たす $l, r$ についてのみ調べる
- $f$ が Monge であるかどうかを判定する
- `output_detail` = `true` のときは、Monge 性を満たさない箇所を具体的に出力する


```cpp
template<class FUNC> bool is_anti_monge
(int N, int M, const FUNC &f, bool upper_triangular = false, bool output_detail = true)
```

- 関数 $f(l, r)$ は、$0 \le l < N$, $0 \le r < M$ に対して定義されるとする
  - `upper_triangular` = `true` のときは、$0 \le l < r \le N$ を満たす $l, r$ についてのみ調べる
- $f$ が anti-Monge であるかどうかを判定する
- `output_detail` = `true` のときは、Monge 性を満たさない箇所を具体的に出力する


## 問題例

- [AtCoder ABC 224 B - Mongeness](https://atcoder.jp/contests/abc224/tasks/abc224_b)

## 参考文献

- tatyam: [Monge の手引き書](https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu)
