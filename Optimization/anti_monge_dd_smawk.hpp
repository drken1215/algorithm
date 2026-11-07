#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "smawk.hpp"

// find shortest path on DAG with anti-Monge cost, by D&D SMAWK, O(NK log N)
// vertex: 0, 1, 2, ..., N
// f(i, j) must be anti-Monge in j >= i + K
// Then, A(j, i) = dp[i] + f(i, j) is also anti-Monge
// slide access NG
template<class VAL, class FUNC> vector<pair<VAL, int>> antiMongeDDSMAWK(int N, const FUNC &f, int K = 1) {
    assert(K >= 1);
    const VAL INF = numeric_limits<VAL>::max() / 4;
    vector<pair<VAL, int>> res(N + 1, make_pair(INF, -1));
    vector<pair<VAL, int>> tmp(N + 1, {INF, -1});
    res[0].first = VAL(0);
    auto relax = [&](int i, int j) -> void {
        if (i != 0 && res[i].second == -1) return;
        VAL val = res[i].first + f(i, j);
        if (val < res[j].first) res[j] = {val, i};
    };
    auto rec = [&](auto &&rec, int left, int right) -> void {
        if (right - left <= 1) return;
        int mid = (left + right) / 2, border = min(right, mid + K - 1);
        rec(rec, left, mid);
        for (int j = mid; j < border; j++) {
            for (int i = left; i < mid && i <= j - K; i++) relax(i, j);
        }
        if (border < right) {
            vector<int> X, Y;
            for (int j = border; j < right; j++) X.push_back(j);
            for (int i = mid - 1; i >= left; i--) {
                if (i == 0 || res[i].second != -1) Y.push_back(i);
            }
            if (!Y.empty()) {
                vector<int> cols(Y.size());
                iota(cols.begin(), cols.end(), 0);
                auto f2 = [&](int j, int k) -> VAL { return res[Y[k]].first + f(Y[k], j); };
                SMAWKRec(X, cols, f2, tmp);
                for (auto x : X) {
                    if (tmp[x].first < res[x].first) res[x] = {tmp[x].first, Y[tmp[x].second]};
                }
            }
        }
        rec(rec, mid, right);
    };
    rec(rec, 0, N + 1);
    return res;
}
