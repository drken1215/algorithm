//
// 木に関する各種クエリ処理
//
// verified:
//   Codeforces Round #614 (Div. 1) C. Xenon's Attack on the Gangs
//     https://codeforces.com/contest/1292/problem/C
//


#include <iostream>
#include <vector>
#include <map>
using namespace std;


using Graph = vector<vector<int>>;
struct RunTree {
    // id[v][w] := the index of node w in G[v]
    vector<map<int,int> > id;

    // num[v][i] := the size of subtree of G[v][i] with parent v
    vector<vector<long long> > num;

    // constructor
    RunTree(const Graph &G, int root = 0) {
        init(G, root);
    }

    // size(u, v) := the size of subtree v with parent u
    long long size(int u, int v) {
        return num[u][id[u][v]];
    }

    // lca(u, v)
    int getLCA(int u, int v) {
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

    // length(u, v)
    long long length(int u, int v) {
        int lca = getLCA(u, v);
        return depth[u] + depth[v] - depth[lca]*2;
    }

    // getParent(v, p) := the parent of v directed for p
    int getParent(int v, int p) {
        if (v == p) return -1;
        int lca = getLCA(v, p);
        if (lca != v) return parent[0][v];
        for (int i = (int)parent.size()-1; i >= 0; --i) {
            if (parent[i][p] != -1 && depth[parent[i][p]] > depth[v]) {
                p = parent[i][p];
            }
        }
        return p;
    }
    
    // rec
    vector<vector<int> > parent;
    vector<int> depth;
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
};



/////////////////////////////////////////////
// solver
/////////////////////////////////////////////

const long long INF = 1LL<<60;
template<class T> inline bool chmax(T& a, T b) { if (a < b) { a = b; return 1; } return 0; }
template<class T> inline bool chmin(T& a, T b) { if (a > b) { a = b; return 1; } return 0; }

using pint = pair<int,int>;
int N;
Graph G;
RunTree rt;

vector<vector<long long> > dp;
long long dprec(int u, int v) {
    if (dp[u][v] != INF) return dp[u][v];
    if (dp[v][u] != INF) return dp[v][u];
    if (u == v) return 0;

    long long res = 0;
    int up = rt.getParent(u, v);
    int vp = rt.getParent(v, u);
    chmax(res, dprec(up, v));
    chmax(res, dprec(vp, u));
    res += rt.size(up, u) * rt.size(vp, v);
    return dp[u][v] = dp[v][u] = res;
}

long long solve() {
    rt.init(G);
    dp.assign(N, vector<long long>(N, INF));
    long long res = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
            chmax(res, dprec(i, j));
        }
    }
    return res;
}

int main() {
    while (scanf("%d", &N) != EOF) {
        G.assign(N, vector<int>());
        for (int i = 0; i < N-1; ++i) {
            int u, v;
            scanf("%d %d", &u, &v);
            --u, --v;
            G[u].push_back(v);
            G[v].push_back(u);
        }
        cout << solve() << endl;
    }
}
