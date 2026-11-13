#define PROBLEM "https://yukicoder.me/problems/no/2764"

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "graph.hpp"
#include "general_li_chao_tree.hpp"

int main() {
    ll N, M, MAX = 1LL << 31, INF = 1LL << 61;
    cin >> N >> M;
    vll W(N);
    cin >> W;
    Graph<ll> G(N);
    REP(i, M) {
        ll u, v, w;
        cin >> u >> v >> w, u--, v--;
        G.add_bidirected_edge(u, v, w);
    }

    auto dijkstra = [&](vll &dp) {
        min_priority_queue<pll> que;
        for (int v = 0; v < N; v++) que.push({dp[v], v});
        while (!que.empty()) {
            auto [cur, v] = que.top(); que.pop();
            if (cur > dp[v]) continue;
            for (auto e : G[v]) {
                if (chmin(dp[e.to], dp[v] + e.val)) {
                    que.push({dp[e.to], e.to});
                }
            }
        }
    };
    vll dp(N, INF), dp2(N, INF);
    dp[0] = 0;
    dijkstra(dp);

    auto make_func = [&](ll a, ll b) {
        return [&, a, b](ll x) -> ll {
            return a * x + b;
        };
    };
    DynamicGeneralLiChaoTree<ll, decltype(make_func(0, 0))> cht(0, MAX);
    for (int v = 0; v < N; v++) cht.insert(make_func(W[v], dp[v]));
    for (int v = 0; v < N; v++) {
        chmin(dp2[v], cht.query(W[v]).first);
    }
    dijkstra(dp2);
    ll res = min(dp[N-1], dp2[N-1]);
    cout << res << endl;
}
