//
// undo つき Union-Find tree
//
// verified:
//   AtCoder ABC 302 Ex - Ball Collector
//     https://atcoder.jp/contests/abc302/tasks/abc302_h
//
//   Codeforces Round 680 (Div.1) C. Team-Building
//     https://codeforces.com/contest/1444/problem/C
//


#include <bits/stdc++.h>
using namespace std;


// Union-Find, we can undo
struct UnionFind {
    // core member
    vector<int> par;
    stack<pair<int,int>> history;
    
    // constructor
    UnionFind() {}
    UnionFind(int n) : par(n, -1) { }
    void init(int n) { par.assign(n, -1); }
    
    // core methods
    int root(int x) {
        if (par[x] < 0) return x;
        else return root(par[x]);
    }
    
    bool same(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y) {
        x = root(x), y = root(y);
        history.emplace(x, par[x]);
        history.emplace(y, par[y]);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y);  // merge technique
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
    
    // debug
    friend ostream& operator << (ostream &s, UnionFind uf) {
        map<int, vector<int>> groups;
        for (int i = 0; i < uf.par.size(); ++i) {
            int r = uf.root(i);
            groups[r].push_back(i);
        }
        for (const auto &it : groups) {
            s << "group: ";
            for (auto v : it.second) s << v << " ";
            s << endl;
        }
        return s;
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void ABC_302_Ex() {
    using pint = pair<int,int>;
    using Graph = vector<vector<int>>;
    
    int N;
    cin >> N;
    vector<int> A(N), B(N);
    for (int i = 0; i < N; ++i) {
        cin >> A[i] >> B[i];
        --A[i], --B[i];
    }
    Graph G(N);
    for (int i = 0; i < N-1; ++i) {
        int u, v;
        cin >> u >> v;
        --u, --v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    
    // Union-Find
    // 外部データの undo もここで実現する
    UnionFind uf(N);
    vector<int> nums(N, 0);  // 各連結成分の辺数
    int cur = 0;  // 現時点での種類数の最大値
    vector<pair<pint,int>> hist;  // 履歴
    vector<int> res(N, 0);  // 答え
    auto insert = [&](int v) -> void {
        int x = uf.root(A[v]), y = uf.root(B[v]);
        hist.push_back(make_pair(pint(x, nums[x]), cur));
        hist.push_back(make_pair(pint(y, nums[y]), cur));
        
    };
    auto erase = [&](int v) -> void {
        
        
        for (int iter = 0; iter < 2; ++iter) {
            nums[hist.top().first.first] = history.top().first.second;
            cur = hist.top().second;
            hist.pop();
        }
    };
    
    // DFS
    auto dfs = [&](auto self, int v, int p) -> void {
        if (uf.same(A))
    };
    dfs(dfs, 0, -1);
    
    for (int v = 1; v < N; ++v) cout << res[v] << " ";
    cout << endl;
}

void CF_680_DIV1_C() {
    using pint = pair<int,int>;
    using Graph = vector<vector<int>>;
    
    cin.tie(0);
    ios::sync_with_stdio(false);

    // 入力
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
    
    // これまでの履歴の削除
    uf.snapshot();

    // すでに二部グラフはダメ
    vector<bool> isbi(K, true);
    long long con = K;
    for (int v = 0; v < N; ++v) {
        if (!isbi[C[v]]) continue;
        if (uf.same(v, v+N)) --con, isbi[C[v]] = false;
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
            if (uf.same(u, u+N) || uf.same(v, v+N)) ok = false;
        }
        if (!ok) --res;

        // 元に戻す
        uf.rollback();
    }
    cout << res << endl;
}


int main() {
    ABC_302_Ex();
    //CF_680_DIV1_C();
}
