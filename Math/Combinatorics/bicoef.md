---
title: 二項係数 (オーソドックスな前処理, $n \le 10^7, r \le 10^7, p \le 10^9$)
documentation_of: ./bicoef.hpp
---

## 概要
二項係数 $\_{n}\mathrm{C}_{r}$ を素数 $p$ で割った余りを求めるライブラリである。素数 $p$ で割った余りの世界を表す構造体をテンプレート引数 `mint` に与えて用いる。

## 参考文献

- けんちょん: [よくやる二項係数 (nCk mod. p)、逆元 (a^-1 mod. p) の求め方](https://drken1215.hatenablog.com/entry/2018/06/08/210000)
