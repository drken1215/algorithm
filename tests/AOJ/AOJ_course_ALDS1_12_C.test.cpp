#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_12_C"

#include <bits/stdc++.h>
using namespace std;

#include "dijkstra.hpp"

int main() {
    int N;
    cin >> N;
    Graph<long long> G(N);
    for (int i = 0; i < N; ++i) {
        int v, num;
        cin >> v >> num;
        for (int j = 0; j < num; j++) {
            int to;
            long long weight;
            cin >> to >> weight;
            G.add_edge(v, to, weight);
        }
    }
    Dijkstra dij(G);
    for (int v = 0; v < N; ++v) cout << v << " " << dij.dp[v] << endl;
}
