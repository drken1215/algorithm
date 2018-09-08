//
// LCA (by Doubling)
//
// verified:
//   AOJ Course GRL_5_C Tree - Lowest Common Ancestor
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_5_C&lang=jp
//


#include <iostream>
#include <vector>
using namespace std;


using Graph = vector<vector<int> >;
struct LCA {
    vector<vector<int> > parent; // parent[d][v] := 2^d-th parent of v
    vector<int> depth;
    LCA() { }
    LCA(const Graph &G, int r = 0) { init(G, r); }
    void init(const Graph &G, int r = 0) {
        int V = (int)G.size();
        int h = 1;
        while ((1<<h) < V) ++h;
        parent.assign(h, vector<int>(V, -1));
        depth.assign(V, -1);
        dfs(G, r, -1, 0);
        for (int i = 0; i+1 < (int)parent.size(); ++i)
            for (int v = 0; v < V; ++v)
                if (parent[i][v] != -1)
                    parent[i+1][v] = parent[i][parent[i][v]];
    }
    void dfs(const Graph &G, int v, int p, int d) {
        parent[0][v] = p;
        depth[v] = d;
        for (auto e : G[v]) if (e != p) dfs(G, e, v, d+1);
    }
    int get(int u, int v) {
        if (depth[u] > depth[v]) swap(u, v);
        for (int i = 0; i < (int)parent.size(); ++i)
            if ( (depth[v] - depth[u]) & (1<<i) )
                v = parent[i][v];
        if (u == v) return u;
        for (int i = (int)parent.size()-1; i >= 0; --i) {
            if (parent[i][u] != parent[i][v]) {
                u = parent[i][u];
                v = parent[i][v];
            }
        }
        return parent[0][u];
    }
};



int main() {
    int N; cin >> N;
    Graph G(N);
    for (int i = 0; i < N; ++i) {
        int num; cin >> num;
        for (int j = 0; j < num; ++j) {
            int c; cin >> c;
            G[i].push_back(c);
            G[c].push_back(i);
        }
    }
    LCA lca(G);
    int Q; cin >> Q;
    for (int q = 0; q < Q; ++q) {
        int u, v; cin >> u >> v;
        cout << lca.get(u, v) << endl;
    }
}
