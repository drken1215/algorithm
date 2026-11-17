#define PROBLEM "https://yukicoder.me/problems/no/2764"

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "graph.hpp"
#include "dijkstra.hpp"
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

    // dp1 step
    Dijkstra dijkstra(G, 0);
    vll dp = dijkstra.get_dp();
    auto make_func = [&](ll a, ll b) {
        return [&, a, b](ll x) -> ll {
            return a * x + b;
        };
    };
    DynamicGeneralLiChaoTree<ll, decltype(make_func(0, 0))> cht(0, MAX);
    for (int v = 0; v < N; v++) {
        cht.insert(make_func(W[v], dp[v]));
    }

    // dp1 -> dp2
    vll dp2(N, INF);
    for (int v = 0; v < N; v++) {
        chmin(dp2[v], cht.query(W[v]).first);
    }

    // dp2 step
    dp2 = dijkstra.solve(G, dp2);
    ll res = min(dp[N-1], dp2[N-1]);
    cout << res << endl;
}
