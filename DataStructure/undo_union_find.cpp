//
// undo つき Union-Find tree
//
// verified:
//   Codeforces Round 680 (Div.1) C. Team-Building
//     https://codeforces.com/contest/1444/problem/C
//


#include <bits/stdc++.h>
using namespace std;


// Union-Find which we can undo
struct UnionFind {
    vector<int> par;
    stack<pair<int,int>> history;
    
    UnionFind() {}
    UnionFind(int n) : par(n, -1) { }
    void init(int n) { par.assign(n, -1); }
    
    int root(int x) {
        if (par[x] < 0) return x;
        else return root(par[x]);
    }
    
    bool issame(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y) {
        x = root(x); y = root(y);
        history.emplace(x, par[x]);
        history.emplace(y, par[y]);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y); // merge technique
        par[x] += par[y];
        par[y] = x;
        return true;
    }
    
    int size(int x) {
        return -par[root(x)];
    }

    // 1-step undo
    void undo() {
        for (int iter = 0; iter < 2; ++iter) {
            par[history.top().first] = history.top().second;
            history.pop();
        }
    }

    // erase history
    void snapshot() {
        while (!history.empty()) history.pop();
    }

    // all rollback
    void rollback() {
        while (!history.empty()) undo();
    }
};

/* debug */
/*
ostream& operator << (ostream& s, UnionFind uf) {
    map<int, vector<int> > res;
    for (int i = 0; i < uf.par.size(); ++i) {
        int r = uf.root(i);
        res[r].push_back(i);
    }
    for (map<int, vector<int> >::iterator it = res.begin(); it != res.end(); ++it) {
        s << endl;
        for (int j = 0; j < (int)it->second.size(); ++j) {
            s << it->second[j] << ", ";
        }
    }
    return s << endl;
}
*/



//////////////////////////////////////////
// solver
//////////////////////////////////////////

using pint = pair<int,int>;
using Graph = vector<vector<int>>;
long long CF680DIV1C(int N, int M, int K, const vector<int> &C, const Graph &G, UnionFind &uf) {
    // これまでの履歴の削除
    uf.snapshot();

    // すでに二部グラフはダメ
    vector<bool> isbi(K, true);
    long long con = K;
    for (int v = 0; v < N; ++v) {
        if (!isbi[C[v]]) continue;
        if (uf.issame(v, v+N)) --con, isbi[C[v]] = false;
    }
    long long res = con * (con - 1) / 2;

    // 辺の色を分類
    map<pint, vector<pint>> ma;
    for (int v = 0; v < N; ++v) {
        for (auto u : G[v]) {
            if (!isbi[C[v]] || !isbi[C[u]] || C[v] == C[u]) continue;
            ma[pint(min(C[v], C[u]), max(C[v], C[u]))].push_back(pint(v, u));
        }
    }

    // 各ペアごとに追加していく
    for (auto it : ma) {
        bool ok = true;
        for (auto e : it.second) {
            int u = e.first, v = e.second;
            uf.merge(u, v+N), uf.merge(u+N, v);
            if (uf.issame(u, u+N) || uf.issame(v, v+N)) ok = false;
        }
        if (!ok) --res;

        // 元に戻す
        uf.rollback();
    }
    return res;
}

int main() {
    cin.tie(0); 
    ios::sync_with_stdio(false);

    int N, M, K;
    cin >> N >> M >> K;
    vector<int> C(N);
    for (int i = 0; i < N; ++i) cin >> C[i], --C[i];
    Graph G(N);
    UnionFind uf(N*2);
    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        --u, --v;
        G[u].push_back(v), G[v].push_back(u);
        if (C[u] == C[v]) uf.merge(u, v+N), uf.merge(u+N, v);
    }
    cout << CF680DIV1C(N, M, K, C, G, uf) << endl;
}
