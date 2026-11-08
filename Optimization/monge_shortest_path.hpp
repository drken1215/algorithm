#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "simple_larsch.hpp"
#include "larsch.hpp"
#include "dd_monotone_minima.hpp"
#include "dd_smawk.hpp"

// Monge Shortest Path Wrapper
// vertex: 0, 1, 2, ..., N
// f(i, j) must be Monge (Then, A(j, i) = dp[i] + f(i, j) is also Monge)
template<class VAL> struct MongeShortestPath {
    VAL INF = numeric_limits<VAL>::max() / 2;
    int CNT_INF = numeric_limits<int>::max() / 2;

    // results
    vector<VAL> dp;
    vector<int> cnt, prev;

    // solver (random access ver)
    template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f
    , const string solver = "simple_larsch") {
        vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
        if (solver == "simple_larsch") {
            SimpleLARSCH<VAL> opt;
            res = opt.solve(N, f);
            dp = opt.dp, cnt = opt.cnt, prev = opt.prev;
        } else if (solver == "larsch") {
            LARSCH<VAL> opt;
            res = opt.solve(N, f);
            dp = opt.dp, cnt = opt.cnt, prev = opt.prev;
        } else if (solver == "monotone") {
            DDMonotoneMinima<VAL> opt;
            res = opt.solve(N, f);
            dp = opt.dp, cnt = opt.cnt, prev = opt.prev;
        } else if (solver == "smawk") {
            DDSMAWK<VAL> opt;
            res = opt.solve(N, f);
            dp = opt.dp, cnt = opt.cnt, prev = opt.prev;
        } else {
            cout << "no solver" << endl;
            assert(true);
        }
        return res;
    }

    // solver (slide access ver)
    template<class STATE, class ADD, class DEL, class GETCOST> vector<pair<VAL, int>> solve
    (int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, const string solver = "simple_larsch") {
        vector<pair<VAL, int>> res;
        if (solver == "simple_larsch") {
            SimpleLARSCH<VAL> opt;
            res = opt.solve(N, ini, add, del, get);
            dp = opt.dp, cnt = opt.cnt, prev = opt.prev;
        } else if (solver == "monotone") {
            DDMonotoneMinima<VAL> opt;
            STATE st = ini;
            int left = 0, right = 0;
            auto f = [&](int i, int j) -> long long {
                while (left < i) del(st, left++);
                while (i < left) add(st, --left);
                while (right < j) add(st, right++);
                while (j < right) del(st, --right);
                return get(st);
            };
            res = opt.solve(N, f);
            dp = opt.dp, cnt = opt.cnt, prev = opt.prev;
        }else {
            cout << "no solver" << endl;
            assert(true);
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
