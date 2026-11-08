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

    // solver (random access ver)
    template<class FUNC> vector<vector<pair<VAL, int>>> solve(int N, const FUNC &f, int D, const string &solver = "smawk") {
        assert(D <= N);
        vector<vector<pair<VAL, int>>> res(N + 1, vector<pair<VAL, int>>(D + 1, make_pair(INF, -1)));
        vector<VAL> dp(N + 1, INF);
        res[0][0].first = 0, dp[0] = 0;
        auto f2 = [&](int i, int j) -> VAL { return (j < i ? dp[j] + f(j, i) : INF); };
        for (int d = 1; d <= D; d++) {
            vector<pair<VAL, int>> tmp;
            if (solver == "smawk") tmp = SMAWK<VAL>(N + 1, N + 1, f2);
            else if (solver == "monotone") tmp = MonotoneMinima<VAL>(N + 1, N + 1, f2);
            else {
                cout << "no solver" << endl;
                assert(true);
            } 
            for (int i = d; i <= N; i++) {
                dp[i] = tmp[i].first;
                res[i][d] = tmp[i];
            }
        }
        return res;
    }

    // solver (slide access ver)
    template<class ADD, class DEL, class GETCOST> vector<vector<pair<VAL, int>>> solve
    (int N, const ADD &add, const DEL &del, const GETCOST &get, int D, const string &solver = "monotone") {
        assert(D <= N);
        vector<vector<pair<VAL, int>>> res(N + 1, vector<pair<VAL, int>>(D + 1, make_pair(INF, -1)));
        vector<VAL> dp(N + 1, INF);
        res[0][0].first = 0, dp[0] = 0;
        int left = 0, right = 0;
        auto f = [&](int i, int j) -> long long {
            while (left < i) del(left++);
            while (i < left) add(--left);
            while (right < j) add(right++);
            while (j < right) del(--right);
            return get();
        };
        auto f2 = [&](int i, int j) -> VAL { return (j < i ? dp[j] + f(j, i) : INF); };
        for (int d = 1; d <= D; d++) {
            vector<pair<VAL, int>> tmp;
            if (solver == "monotone") tmp = MonotoneMinima<VAL>(N + 1, N + 1, f2);
            else {
                cout << "no solver" << endl;
                assert(true);
            } 
            for (int i = d; i <= N; i++) {
                dp[i] = tmp[i].first;
                res[i][d] = tmp[i];
            }
        }
        return res;
    }
};
