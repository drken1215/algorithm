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
template<class VAL> struct DDSMAWK {
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
        vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
        vector<pair<VAL, int>> tmp(N + 1);
        dp[0] = 0, cnt[0] = 0, res[0].first = VAL(0);
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
