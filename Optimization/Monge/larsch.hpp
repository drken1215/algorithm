#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// LARSCH, in O(N)
// vertex: 0, 1, 2, ..., N
// f(i, j) must be Monge (Then, A(j, i) = dp[i] + f(i, j) is also Monge)
template<class VAL, class FUNC> struct SubLARSCH {
    struct RowReduce;
    struct ColReduce;
    struct ColMap {
        const ColMap* parent = nullptr;
        const vector<int> *v = nullptr;
        constexpr int mapping(int i) const {
            int x = (v ? (*v)[i] : i);
            return (parent ? parent->mapping(x) : x);
        }
    };
    struct Eval {
        const FUNC *f = nullptr;
        long long a = 1, b = 0;  // row = a * i + b
        const ColMap *cm = nullptr;
        constexpr VAL operator () (int i, int j) const {
            int i2 = int(a * i + b), j2 = (cm ? cm -> mapping(j) : j);
            return (*f)(i2, j2);
        }
    };
    struct RowReduce {
        int N;
        Eval e;
        int cur_row = 0, state = 0;
        unique_ptr<ColReduce> rec;
        constexpr RowReduce(int N_, const Eval &e_) : N(N_), e(e_) {
            int M = N / 2;
            if (M > 0) {
                Eval e2 = e;
                e2.a = e.a * 2, e2.b = e.a + e.b;
                rec = make_unique<ColReduce>(M, e2);
            }
        }
        constexpr void reset() {
            cur_row = 0, state = 0;
            if (rec) rec->reset();
        }
        int get_argmin() {
            int i = cur_row++;
            if (!(i & 1)) {
                int prev = state, next = (i + 1 == N ? N - 1 : rec->get_argmin());
                state = next;
                int res = prev;
                for (int j = prev + 1; j <= next; j++) {
                    if (e(i, res) > e(i, j)) res = j;
                }
                return res;
            } else {
                return (e(i, state) <= e(i, i) ? state : i);
            }
        }
    };
    struct ColReduce {
        int N;
        Eval e;
        int cur_row = 0;
        vector<int> cols;
        ColMap cm;
        RowReduce rec;
        constexpr ColReduce(int N_, const Eval &e_) 
        : N(N_), e(e_), cols(), cm{e.cm, &cols}, rec(N_, Eval{e.f, e.a, e.b, &cm}) {
            cols.reserve(N);
        } 
        constexpr void reset() {
            cur_row = 0;
            cols.clear();
            rec.reset();
        }
        constexpr void push_col(int i, int j) {
            while (!cols.empty()) {
                int siz = (int)cols.size();
                if (siz == i) break;
                int last = cols.back();
                if (e(siz - 1, last) > e(siz - 1, j)) cols.pop_back();
                else break;
            }
            if ((int)cols.size() != N) cols.emplace_back(j);
        }
        constexpr int get_argmin() {
            int i = cur_row++;
            if (i == 0) {
                cols.clear();
                cols.emplace_back(0);
            } else {
                push_col(i, i * 2 - 1);
                push_col(i, i * 2);
            }
            return cols[rec.get_argmin()];
        }
    };

    FUNC f;
    ColMap root_cm;
    Eval root_eval;
    unique_ptr<RowReduce> base;
    explicit SubLARSCH(int N_, const FUNC &f_) 
    : f(std::move(f_)), root_cm{nullptr, nullptr}, root_eval{&f, 1, 0, &root_cm} {
        base = make_unique<RowReduce>(N_, root_eval);
    }
    constexpr void reset() {
        base->reset();
    }
    constexpr int get_argmin() {
        return base->get_argmin();
    }
};

template<class VAL> struct LARSCH {
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
        dp[0] = 0, cnt[0] = 0, res[0].first = 0;
        auto f2 = [&](int i, int j) -> VAL {
            i++;
            if (i <= j) return INF;
            else return dp[j] + f(j, i);
        };
        SubLARSCH<VAL, decltype(f2)> lar(N, f2);
        for (int r = 1; r <= N; r++) {
            int l = lar.get_argmin();
            dp[r] = dp[l] + f(l, r);
            prev[r] = l;
            cnt[r] = cnt[l] + 1;
        }
        for (int r = 1; r <= N; r++) res[r] = {dp[r], prev[r]};
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
