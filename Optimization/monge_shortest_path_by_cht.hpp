#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "general_li_chao_tree.hpp"

// find shortest path on DAG with Pairwise Single Crossing cost, by Li Chao Tree, O(N log N)
// vertex: 0, 1, 2, ..., N
// A(j, i) = dp[i] + f(i, j) must be Pairwise Single Crossing (must be in all 0 <= i, j <= N, including i >= j)
// Attention: Even if f is Pairwise Single Crossing, A is not necessarily Pairwise Single Crossing
// edges: (i, j), j >= start[i] (i + 1 <= start[i] <= N + 1), start is non-decreasing
template<class VAL, class FUNC> vector<pair<VAL, int>> MongeShortestPathByCHT(int N, const FUNC &f, const vector<int> &start) {
    assert((int)start.size() == N);
    for (int i = 0; i < N; i++) {
        assert(i + 1 <= start[i] && start[i] <= N + 1);
        if (i > 0) assert(start[i - 1] <= start[i]);
    }
    VAL INF = numeric_limits<VAL>::max() / 2;
    vector<pair<VAL, int>> dp(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
    auto make_func = [&](int i) {
        return [&, i](long long j) -> VAL { return dp[i].first + f(i, j); };
    };
    GeneralLiChaoTree<VAL, decltype(make_func(0))> cht(1, N + 1, INF);
    dp[0].first = 0;
    int i = 0;
    for (int j = 1; j <= N; j++) {
        while (i < N && j >= start[i]) {
            cht.insert(make_func(i));
            i++;
        }
        dp[j] = cht.query(j);
    };
    return dp;
}
template<class VAL, class FUNC> vector<pair<VAL, int>> MongeShortestPathByCHT(int N, const FUNC &f) {
    vector<int> start(N);
    for (int i = 0; i < N; i++) start[i] = i + 1;
    return MongeShortestPathByCHT<VAL>(N, f, start);
}
