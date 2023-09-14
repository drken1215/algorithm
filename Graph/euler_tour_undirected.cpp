//
// 無向オイラー路を求める
//
// verified:
//   Codeforces 554 DIV2 E. Neko and Flashback
//     https://codeforces.com/contest/1152/problem/E
//


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


// edge class (for network-flow)
template<class FLOWTYPE> struct Edge {
    int rev, from, to;
    FLOWTYPE cap, icap;
    Edge(int r, int f, int t, FLOWTYPE c) : rev(r), from(f), to(t), cap(c), icap(c) {}
    friend ostream& operator << (ostream& s, const Edge& E) {
        if (E.cap > 0) return s << E.from << "->" << E.to << '(' << E.cap << ')';
        else return s;
    }
};

// graph class (for network-flow)
template<class FLOWTYPE> struct Graph {
    vector<vector<Edge<FLOWTYPE> > > list;
    
    Graph(int n = 0) : list(n) { }
    void init(int n = 0) { list.clear(); list.resize(n); }
    void reset() { for (int i = 0; i < (int)list.size(); ++i) for (int j = 0; j < list[i].size(); ++j) list[i][j].cap = list[i][j].icap; }
    inline vector<Edge<FLOWTYPE> >& operator [] (int i) { return list[i]; }
    inline const size_t size() const { return list.size(); }
    
    inline Edge<FLOWTYPE> &redge(Edge<FLOWTYPE> e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    
    void addedge(int from, int to, FLOWTYPE cap) {
        list[from].push_back(Edge<FLOWTYPE>((int)list[to].size(), from, to, cap));
        list[to].push_back(Edge<FLOWTYPE>((int)list[from].size() - 1, to, from, 0));
    }
    
    void add_undirected_edge(int from, int to, FLOWTYPE cap) {
        list[from].push_back(Edge<FLOWTYPE>((int)list[to].size(), from, to, cap));
        list[to].push_back(Edge<FLOWTYPE>((int)list[from].size() - 1, to, from, cap));
    }

    /*
    // debug
    friend ostream& operator << (ostream& s, const Graph& G) {
        s << endl; for (int i = 0; i < G.size(); ++i) { s << i << " : " << G.list[i] << endl; }return s;
    }
    */
};

struct UnionFind {
    vector<int> par;
    
    UnionFind(int n) : par(n, -1) { }
    void init(int n) { par.assign(n, -1); }
    
    int root(int x) {
        if (par[x] < 0) return x;
        else return par[x] = root(par[x]);
    }
    
    bool issame(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y) {
        x = root(x); y = root(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y); // merge technique
        par[x] += par[y];
        par[y] = x;
        return true;
    }
    
    int size(int x) {
        return -par[root(x)];
    }
};

struct Euler {
    vector<int> nextv, walk;
    Euler() {}
    void dfs(Graph<int> &G, int v) {
        for (int &i = nextv[v]; i < (int)G[v].size(); ++i) {
            auto &e = G[v][i];
            if (e.cap == 0) continue;
            --e.cap;
            --G.redge(e).cap;
            dfs(G, e.to);
        }
        walk.push_back(v);
    }
    vector<vector<int> > solve(Graph<int> &G) {
        vector<vector<int> > res;
        int V = (int)G.size();
        UnionFind uf(V);
        for (int v = 0; v < V; ++v) for (auto e : G[v]) uf.merge(e.from, e.to);
        vector<vector<int> > comps(V);
        for (int v = 0; v < V; ++v) {
            int r = uf.root(v);
            comps[r].push_back(v);
        }
        nextv.assign(V, 0);
        for (int v = 0; v < V; ++v) {
            if (comps[v].empty()) continue;
            vector<int> odds;
            for (auto v2 : comps[v]) if (G[v2].size() & 1) odds.push_back(v2);
            if (odds.size() > 2) return vector<vector<int> >();
            walk.clear();
            int start = (!odds.empty() ? odds[0] : v);
            dfs(G, start);
            reverse(walk.begin(), walk.end());
            res.push_back(walk);
        }
        return res;
    }
};


void solve(const vector<int> &a, const vector<int> &b) {
    vector<int> vals;
    for (int i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) {
            cout << -1 << endl;
            return;
        }
        vals.push_back(a[i]);
        vals.push_back(b[i]);
    }
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());

    Graph<int> G((int)vals.size());
    for (int i = 0; i < (int)a.size(); ++i) {
        int ai = lower_bound(vals.begin(), vals.end(), a[i]) - vals.begin();
        int bi = lower_bound(vals.begin(), vals.end(), b[i]) - vals.begin();
        G.add_undirected_edge(ai, bi, 1);
    }
    Euler el;
    auto res = el.solve(G);
    if (res.empty() || res.size() > 1) {
        cout << -1 << endl;
        return;
    }
    for (int i = 0; i < (int)res[0].size(); ++i) {
        if (i) cout << " ";
        cout << vals[res[0][i]];
    }
    cout << endl;
}

int main() {
    int N;
    while (cin >> N) {
        vector<int> a(N-1), b(N-1);
        for (int i = 0; i < N-1; ++i) cin >> a[i];
        for (int i = 0; i < N-1; ++i) cin >> b[i];
        solve(a, b);
    }
}
