//
// 木を走査して、さまざまな情報を求める
//
// verified:
//   Codeforces Round #614 (Div. 1) C. Xenon's Attack on the Gangs
//     https://codeforces.com/contest/1292/problem/C
//
//   CodeQUEEN 決勝 C - Path Intersection
//     https://atcoder.jp/contests/codequeen2023-final-open/tasks/codequeen2023_final_c
//


#include <bits/stdc++.h>
using namespace std;


// Run Tree
using Graph = vector<vector<int>>;
struct RunTree {
    // id[v][w] := the index of node w in G[v]
    vector<map<int, int>> id;

    // num[v][i] := the size of subtree of G[v][i] with parent v
    vector<vector<long long>> num;
    
    // for finding lca
    vector<vector<int>> parent;
    vector<int> depth;

    // constructor
    RunTree() {}
    RunTree(const Graph &G, int root = 0) {
        init(G, root);
    }
    
    // init
    void init(const Graph &G, int root = 0) {
        int N = (int)G.size();
        id.assign(N, map<int,int>());
        num.assign(N, vector<long long>());
        for (int v = 0; v < N; ++v) num[v].assign((int)G[v].size(), 0);
        int V = (int)G.size();
        int h = 1;
        while ((1<<h) < N) ++h;
        parent.assign(h, vector<int>(N, -1));
        depth.assign(N, -1);
        rec(G, root);
        for (int i = 0; i+1 < (int)parent.size(); ++i)
            for (int v = 0; v < V; ++v)
                if (parent[i][v] != -1)
                    parent[i+1][v] = parent[i][parent[i][v]];
    }

    // size(u, v) := the size of subtree v with parent u
    long long size(int u, int v) {
        return num[u][id[u][v]];
    }

    // lca(u, v)
    int get_lca(int u, int v) {
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

    // dist(u, v)
    long long get_dist(int u, int v) {
        int lca = get_lca(u, v);
        return depth[u] + depth[v] - depth[lca]*2;
    }

    // get_parent(v, p) := the parent of v directed for p
    int get_parent(int v, int p) {
        if (v == p) return -1;
        int lca = get_lca(v, p);
        if (lca != v) return parent[0][v];
        for (int i = (int)parent.size()-1; i >= 0; --i) {
            if (parent[i][p] != -1 && depth[parent[i][p]] > depth[v]) {
                p = parent[i][p];
            }
        }
        return p;
    }
    
    // rec
    int rec(const Graph &G, int v, int p = -1, int d = 0) {
        int p_index = -1;
        int sum = 1;
        parent[0][v] = p;
        depth[v] = d;
        for (int i = 0; i < (int)G[v].size(); ++i) {
            int ch = G[v][i];
            id[v][ch] = i;
            if (ch == p) {
                p_index = i;
                continue;
            }
            int s = rec(G, ch, v, d+1);
            num[v][i] = s;
            sum += s;
        }
        if (p_index != -1) num[v][p_index] = (int)G.size() - sum;
        return sum;
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void Codeforces_614_C() {
    int N;
    scanf("%d", &N);
    Graph G(N);
    for (int i = 0; i < N-1; ++i) {
        int u, v;
        scanf("%d %d", &u, &v);
        --u, --v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    
    // Run Tree
    RunTree rt(G);
    vector<vector<long long>> dp(N, vector<long long>(N, -1));
    
    // メモ化再帰
    auto rec = [&](auto self, int u, int v) -> long long {
        if (dp[u][v] != -1) return dp[u][v];
        if (dp[v][u] != -1) return dp[v][u];
        if (u == v) return 0;

        long long res = 0;
        int up = rt.get_parent(u, v), vp = rt.get_parent(v, u);
        res = max(res, self(self, up, v));
        res = max(res, self(self, vp, u));
        res += rt.size(up, u) * rt.size(vp, v);
        return dp[u][v] = dp[v][u] = res;
    };
    
    // 集計
    long long res = 0;
    for (int i = 0; i < N; ++i) for (int j = i+1; j < N; ++j) {
        res = max(res, rec(rec, i, j));
    }
    cout << res << endl;
}

void CodeQUEEN_D() {
    int N, S, T;
    cin >> N >> S >> T;
    --S, --T;
    Graph G(N);
    for (int i = 0; i < N-1; ++i) {
        int u, v;
        cin >> u >> v;
        --u, --v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    
    RunTree rt(G, S);
    for (int v = 0; v < N; ++v) {
        int lca = rt.get_lca(v, T);
        cout << rt.get_dist(lca, v) + 1 << endl;
    }
}


int main() {
    Codeforces_614_C();
    //CodeQUEEN_D();
}

