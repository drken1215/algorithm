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
template<class VAL> struct DDMonotoneMinima {
    VAL INF = numeric_limits<VAL>::max() / 2;
    int CNT_INF = numeric_limits<int>::max() / 2;

    // results
    vector<VAL> dp;
    vector<int> cnt, prev;

    // solver
    template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f) {
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;
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
        for (int r = 1; r <= N; r++) {
            dp[r] = res[r].first;
            prev[r] = res[r].second;
            if (prev[r] != -1) cnt[r] = cnt[prev[r]] + 1;
        }
        return res;
    }

    // reconstruct the shortest path
    vector<int> reconstruct() {
        int N = (int)dp.size() - 1;
        vector<int> path;
        for (int v = N; v > 0; v = prev[v]) path.emplace_back(v);
        path.emplace_back(0);
        reverse(path.begin(), path.end());
        return path;
    }
};
