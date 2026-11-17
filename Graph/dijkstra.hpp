#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

#include "graph.hpp"

// Dijkstra 法
template<class GRAPH> struct Dijkstra {
    using VAL = typename GRAPH::value_type;
    using EDGE = typename GRAPH::edge_type;
    VAL INF = numeric_limits<VAL>::max() / 2;

    // results
    vector<VAL> dp;  // 最短路長
    vector<int> prev;  // 最短路の直前の頂点
    GRAPH dag;  // 最短路の復元構造を表す DAG

    // constructors and getters
    Dijkstra(const GRAPH &G, int s = 0) {
        solve(G, s);
    }
    vector<VAL> get_dp() const { return dp; }
    vector<VAL> get_dist() const { return dp; }
    vector<int> get_prev() const { return prev; }
    bool is_reachable(int v) {
        assert(v >= 0 && v < (int)dp.size());
        return !dag[v].empty();
    }

    // solver
    vector<VAL> solve(const GRAPH &G, int s) {
        // from s
        dp.assign((int)G.size(), INF);
        prev.assign((int)G.size(), -1);
        dag.assign((int)G.size());
        dp[s] = 0;
        priority_queue<pair<VAL, int>, vector<pair<VAL, int>>, greater<pair<VAL, int>>> que;
        que.push({0, s});
        return inner_solve(G, que);
    }
    vector<VAL> solve(const GRAPH &G, const vector<VAL> &ini) {
        // from initial dp table ini
        dp.assign((int)G.size(), INF);
        prev.assign((int)G.size(), -1);
        priority_queue<pair<VAL, int>, vector<pair<VAL, int>>, greater<pair<VAL, int>>> que;
        for (int v = 0; v < (int)G.size(); v++) {
            assert(ini[v] >= 0);
            if (ini[v] < INF) {
                dp[v] = ini[v];
                que.push({dp[v], v});
            }
        }
        return inner_solve(G, que);
    }
    vector<VAL> inner_solve(const GRAPH &G, priority_queue<pair<VAL, int>, vector<pair<VAL, int>>, greater<pair<VAL, int>>> &que) {
        while (!que.empty()) {
            auto [cur, v] = que.top();
            que.pop();
            if (cur > dp[v]) continue;
            for (auto e : G[v]) {
                if (dp[e.to] > dp[v] + e.val) {
                    dp[e.to] = dp[v] + e.val;
                    prev[e.to] = v;
                    dag[e.to].clear();
                    dag.add_edge(e.to, v, e.val);
                    que.push({dp[e.to], e.to});
                } else if (dp[e.to] == dp[v] + e.val) {
                    dag.add_edge(e.to, v, e.val);
                }
            }
        }
        return dp;
    }

    // reconstruct the shortest path
    vector<EDGE> reconstruct(int s, int t) const {
        vector<EDGE> res;
        int v = t;
        while (!dag[v].empty() && v != s) {
            EDGE e = dag[v][0];
            res.push_back(EDGE{e.to, e.from, e.val});
            v = e.to;
        }
        reverse(res.begin(), res.end());
        return res;
    }
};
