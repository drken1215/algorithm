#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "monotone_minima.hpp"
#include "smawk.hpp"

// find the d-edges shortest path for d = 1, 2, ..., D, in O(ND)
// vertex: 0, 1, 2, ..., N
// f(i, j) must be Monge (Then, A(j, i) = dp[i] + f(i, j) is also Monge)
// you can choose SMAWK (random access) or Monotone Minima (slide access)
template<class VAL> struct MongeShortestPathWithDEdges {
    VAL INF = numeric_limits<VAL>::max() / 2;

    // solver
    template<class FUNC> vector<VAL> solve(int N, const FUNC &f, int D, const string &solver = "smawk") {
        assert(D <= N);
        vector<VAL> res{VAL(0)}, dp(N + 1, INF);
        dp[0] = 0;
        auto f2 = [&](int i, int j) -> VAL { return (j < i ? dp[j] + f(j, i) : INF); };
        for (int d = 1; d <= D; d++) {
            vector<pair<VAL, int>> tmp;
            if (solver == "smawk") tmp = SMAWK<VAL>(N + 1, N + 1, f2);
            else if (solver == "monotone") tmp = MonotoneMinima<VAL>(N + 1, N + 1, f2);
            for (int i = d; i <= N; i++) dp[i] = tmp[i].first;
            res.emplace_back(dp[N]);
        }
        return res;
    }
};
