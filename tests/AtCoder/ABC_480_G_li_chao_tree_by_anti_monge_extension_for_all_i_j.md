---
title: AtCoder ABC 480 G - Zonal Score Minimization (全域 anti-Monge 拡張方法)
documentation_of: ./ABC_480_G_li_chao_tree_by_anti_monge_extension_for_all_i_j.cpp
---

## 問題概要
長さ $N$ の数列 $A_{0}, A_{1}, \dots, A_{N-1}$ をいくつかの区間に分割する。各区間の長さは $2$ 以上でなければならない。

また、各区間のスコアを「区間内部の値の最大値と最小値の差」として定める。

各区間のスコアの総和として考えられる最小値を求めよ。


### 制約

- $2 \le N \le 2 \times 10^{5}$


## 解法メモ
$0 \le i, j \le N-1$ かつ $j \ge i + 2$ を満たす $(i, j)$ に対して、次の関数 $f$ を定める。

----
$f(i, j)$ := $A_{i}, A_{i+1}, \dots, A\_{j-1}$ の最大値と最小値の差

----

このとき、次の DP によって最適解を求めることができる。

`dp[j] = min_{i <= j-2} (dp[i] + f(i, j))`

ここで、$0 \le i, j \le N-1$ かつ $j \ge i + 2$ を満たす領域においては、$f(i, j)$ は **anti-Monge** である。このことを活かして以下の解法が考えられる。

- General Li Chao Tree を用いる解法
  - $j \ge i + 2$ を満たす $i, j$ だけでなく、$0 \le i, j \le N-1$ を満たすすべての $i, j$ に対して anti-Monge となるように $f$ を拡張する解法
  - Li Chao Tree に関数 $f_i$ を挿入するときに、定義域を $[i+2, N]$ に限定して区間挿入する解法
- 分割統治 SMAWK 法を用いる解法

本コードは、このうちの「すべての $i, j$ に対して anti-Monge となるように $f$ を拡張する解法」である。具体的には、十分大きな定数 $M$ を用いて、次のコードのようにすればよい。

```cpp
// Li Chao Tree は、通常、全定義域に anti-Monge 性を拡張する必要がある！！
auto func = [&](int l, int r) -> long long { 
    if (r < l + 2) return -M * ((l + 1 - r) * 2 + 1);
    else return stma.get(l, r) - stmi.get(l, r);
};
```

