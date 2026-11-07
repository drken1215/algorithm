#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "monotone_minima.hpp"

// find shortest path on DAG by D&D Monotone Minima, O(N (log N)^2)
// vertex: 0, 1, 2, ..., N
// A(j, i) = dp[i] + f(i, j) must be monotone (argmin is not decreasing)
// Attention: Even if f is monotone, A is not necessarily monotone
// slide access OK
template<class VAL, class FUNC> vector<pair<VAL, int>> DDMonotoneMinima(int N, const FUNC &f) {
    vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
    res[0].first = VAL(0);
    auto f2 = [&](int i, int j) -> VAL { return res[j].first + f(j, i); };
    auto rec = [&](auto &&rec, int left, int right) -> void {
        if (right - left <= 1) return;
        int mid = (left + right) / 2;
        rec(rec, left, mid);
        MonotoneMinimaRec(mid, right, left, mid, f2, res);
        rec(rec, mid, right);
    };
    rec(rec, 0, N + 1);
    return res;
}
