#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "smawk.hpp"

// find shortest path on DAG by D&D SMAWK, O(N log N)
// vertex: 0, 1, 2, ..., N
// A(j, i) = dp[i] + f(i, j) must be totally monotone
// Attention: Even if f is totally monotone, A is not necessarily totally monotone
// slide access NG
template<class VAL, class FUNC> vector<pair<VAL, int>> DDSMAWK(int N, const FUNC &f) {
    vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
    vector<pair<VAL, int>> tmp(N + 1);
    res[0].first = VAL(0);
    auto f2 = [&](int i, int j) -> VAL { return res[j].first + f(j, i); };
    auto rec = [&](auto &&rec, int left, int right) -> void {
        if (right - left <= 1) return;
        int mid = (left + right) / 2;
        vector<int> X(right - mid), Y(mid - left);
        for (int i = mid; i < right; i++) X[i - mid] = i;
        for (int j = left; j < mid; j++) Y[j - left] = j;
        rec(rec, left, mid);
        SMAWKRec(X, Y, f2, tmp);
        for (auto x : X) if (tmp[x].first < res[x].first) res[x] = tmp[x];
        rec(rec, mid, right);
    };
    rec(rec, 0, N + 1);
    return res;
}
