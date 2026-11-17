#define PROBLEM "https://judge.yosupo.jp/problem/shortest_path"

#include <bits/stdc++.h>
using namespace std;

#include "dijkstra.hpp"

int main() {
    int N, M, s, t;
    cin >> N >> M >> s >> t;
    Graph<long long> G(N);
    for (int i = 0; i < M; ++i) {
        int a, b, c;
        cin >> a >> b >> c;
        G.add_edge(a, b, c);
    }
    Dijkstra dij(G, s);
    const auto &path = dij.reconstruct(s, t);
    
    if (!dij.is_reachable(t)) cout << -1 << endl;
    else {
        cout << dij.dp[t] << " " << path.size() << endl;
        for (auto e : path) {
            cout << e.from << " " << e.to << endl;
        }
    }
}
