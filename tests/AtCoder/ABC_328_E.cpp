#include <bits/stdc++.h>
using namespace std;

#include "union_find.hpp"
#include "next_combination.hpp"


int main() {
    struct Edge {
        int from, to;
        long long weight;
        Edge(int f = 0, int t = 0, long long w = 0) : from(f), to(t), weight(w) {}
    };
    long long N, M, K, res = 1LL << 60;
    cin >> N >> M >> K;
    vector<Edge> edges(M);
    for (int i = 0; i < M; ++i) {
        long long u, v, w;
        cin >> u >> v >> w;
        --u, --v;
        edges[i] = Edge(u, v, w);
    }
    long long bit = (1LL << (N-1)) - 1;
    do {
        long long tmp = 0;
        bool ok = true;
        UnionFind uf(N);
        for (int i = 0; i < M; ++i) {
            if (bit >> i & 1) {
                auto [f, t, w] = edges[i];
                if (uf.same(f, t)) ok = false;
                uf.merge(f, t);
                tmp += w;
            }
        }
        if (ok) res = min(res, tmp % K);
    } while (next_combination(bit, M));
    cout << res << endl;
}
