#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "monge_shortest_path.hpp"
#include "alien_trick.hpp"

// Alien DP
// find the d-edges shortest path from vertex 0 on DAG by Lagrange relaxation
// vertex: 0, 1, 2, ..., N, f(i, j) must be Monge
template<class VAL> struct AlienDP {
    // inner values
    MongeShortestPath<VAL> msp;
    AlienTrick<VAL> at;

    // solver (random access ver.)
    template<class FUNC> pair<VAL, VAL> solve(int N, const FUNC &f, int D, const string &solver = "simple_larsch") {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            vector<pair<VAL, int>> obj = msp.solve(N, fg, solver);
            auto cnt = msp.cnt[N];
            auto res = obj.back().first - lambda * cnt;
            return {res, cnt};
        };
        return at.solve(lag, D);
    }

    // solver (slide access ver.)
    template<class STATE, class ADD, class DEL, class GETCOST> pair<VAL, int> solve
    (int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, int D, const string &solver = "simple_larsch") {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto get_lag = [&](STATE &st) -> VAL { return get(st) + lambda; };
            vector<pair<VAL, int>> obj = msp.solve(N, ini, add, del, get_lag, solver);
            auto cnt = msp.cnt[N];
            auto res = obj.back().first - lambda * cnt;
            return {res, cnt};
        };
        return at.solve(lag, D);
    }

    // debugger
    template<class FUNC> void debug(int N, const FUNC &f, VAL start = 0, VAL goal = 10) {
        auto lag = [&](VAL lambda) -> void {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            const auto &dp = msp.solve(N, fg);
            auto cnt = msp.cnt[N];
            auto res = dp.back().first - lambda * cnt;
            cout << lambda << ": " << res << "(" << cnt << "); (";
            for (int i = 0; i < (int)dp.size(); i++) {
                if (i) cout << " ";
                cout << "<" << dp[i].first << ", " << dp[i].second << ">";
            }
            cout << ")" << endl;
        };
        for (VAL lambda = start; lambda <= goal; lambda++) lag(lambda);
    }
};
