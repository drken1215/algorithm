#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// noshi's simplified LARSCH
// find shortest path from vertex 0 on DAG in O(N log N)
// vertex: 0, 1, 2, ..., N
// f(i, j) must be Monge (Then, A(j, i) = dp[i] + f(i, j) is also Monge)
// slide access OK
template<class VAL> struct SimpleLARSCH {
    VAL INF = numeric_limits<VAL>::max() / 2;
    int CNT_INF = numeric_limits<int>::max() / 2;

    // results
    vector<VAL> dp;
    vector<int> cnt, prev;

    // solver (random access ver)
    template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f, bool minimize_cnt = true) {
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;

        auto relax = [&](int l, int r) -> void {
            VAL val = dp[l] + f(l, r);
            int cn = cnt[l] + 1;
            if (dp[r] > val || (dp[r] == val && (minimize_cnt ? cn < cnt[r] : cn > cnt[r]))) {
                dp[r] = val;
                cnt[r] = cn;
                prev[r] = l;
            }
        };
        auto rec = [&](auto &&rec, int l, int r) -> void {
            if (r - l <= 1) return;
            int m = (l + r) / 2;
            for (int k = prev[l]; k <= prev[r]; k++) relax(k, m);
            rec(rec, l, m);
            for (int k = l + 1; k <= m; k++) relax(k, r);
            rec(rec, m, r);
        };

        if (N > 0) relax(0, N), rec(rec, 0, N);
        vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
        res[0].first = VAL(0);
        for (int i = 1; i <= N; i++) res[i] = {dp[i], prev[i]};
        return res;
    }

    // solver (slide access ver)
    template<class STATE, class ADD, class DEL, class GETCOST> vector<pair<VAL, int>> solve
    (int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, bool minimize_cnt = true) {
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;

        STATE win[2] = { ini, ini };
        int cl[2] = {0, 0}, cr[2] = {0, 0};

        auto move_cursor = [&](int c, int l, int r) {
            while (cr[c] < r) add(win[c], cr[c]++);
            while (cl[c] > l) add(win[c], --cl[c]);
            while (cr[c] > r) del(win[c], --cr[c]);
            while (cl[c] < l) del(win[c], cl[c]++);
        };
        auto relax = [&](int c, int l, int r) {
            move_cursor(c, l, r);
            VAL val = dp[l] + get(win[c]);
            int cn = cnt[l] + 1;
            if (dp[r] > val || (dp[r] == val && (minimize_cnt ? cn < cnt[r] : cn > cnt[r]))) {
                dp[r] = val;
                cnt[r] = cn;
                prev[r] = l;
            }
        };
        auto rec = [&](auto &&rec, int l, int r) -> void {
            if (r - l <= 1) return;
            int m = (l + r) / 2;
            for (int k = prev[l]; k <= prev[r]; k++) relax(0, k, m);
            rec(rec, l, m);
            for (int k = l + 1; k <= m; k++) relax(1, k, r);
            rec(rec, m, r);
        };

        if (N > 0) relax(0, 0, N), rec(rec, 0, N);
        vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
        res[0].first = VAL(0);
        for (int i = 1; i <= N; i++) res[i] = {dp[i], prev[i]};
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
