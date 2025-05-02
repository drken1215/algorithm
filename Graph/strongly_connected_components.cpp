//
// 強連結成分分解
//
// verified:
//   Yosupo Judge - Strongly Connected Components
//     https://judge.yosupo.jp/problem/scc
//
//   AOJ Course GRL_3_C - 強連結成分分解
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_3_C&lang=ja
//


#include <bits/stdc++.h>
using namespace std;


// Edge Class
template<class T> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1), val(-1) { }
    Edge(int f, int t, T v = -1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& E) {
        return s << E.from << "->" << E.to;
    }
};

// graph class
template<class T> struct Graph {
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    
    Graph(int n = 0) : list(n), reversed_list(n) { }
    void init(int n = 0) {
        list.assign(n, vector<Edge<T>>());
        reversed_list.assign(n, vector<Edge<T>>());
    }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    const vector<Edge<T>> &get_rev_edges(int i) const { return reversed_list[i]; }
    const size_t size() const { return list.size(); }
        
    void add_edge(int from, int to, T val = -1) {
        list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }
    
    void add_bidirected_edge(int from, int to, T val = -1) {
        list[from].push_back(Edge(from, to, val));
        list[to].push_back(Edge(to, from, val));
        reversed_list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }

    friend ostream &operator << (ostream &s, const Graph &G) {
        s << endl;
        for (int i = 0; i < G.size(); ++i) {
            s << i << " -> ";
            for (const auto &e : G[i]) s << e.to << " ";
            s << endl;
        }
        return s;
    }
};

// strongly connected components decomposition
template<class T> struct SCC {
    // results
    vector<int> cmp;
    vector<vector<int>> groups;
    Graph<T> dag;
    
    // intermediate results
    vector<bool> seen;
    vector<int> vs, rvs;
    
    // constructor
    SCC() { }
    SCC(const Graph<T> &G) { 
        solve(G);
    }
    void init(const Graph<T> &G) { 
        solve(G);
    }

    // getter, compressed dag（v: node-id of compressed dag)
    int get_size(int v) const {
        return groups[v].size();
    }
    vector<int> get_group(int v) const {
        return groups[v];
    }

    // solver
    void dfs(const Graph<T> &G, int v) {
        seen[v] = true;
        for (const auto &e : G[v]) if (!seen[e.to]) dfs(G, e.to);
        vs.push_back(v);
    }
    void rdfs(const Graph<T> &G, int v, int k) {
        seen[v] = true;
        cmp[v] = k;
        for (const auto &e : G.get_rev_edges(v)) if (!seen[e.to]) rdfs(G, e.to, k);
        rvs.push_back(v);
    }
    void reconstruct(const Graph<T> &G) {
        dag.init((int)groups.size());
        set<pair<int,int>> new_edges;
        for (int i = 0; i < (int)G.size(); ++i) {
            int u = cmp[i];
            for (const auto &e : G[i]) {
                int v = cmp[e.to];
                if (u == v) continue;
                if (!new_edges.count({u, v})) {
                    dag.add_edge(u, v);
                    new_edges.insert({u, v});
                }
            }
        }
    }
    void solve(const Graph<T> &G) {
        // first dfs
        seen.assign((int)G.size(), false);
        vs.clear();
        for (int v = 0; v < (int)G.size(); ++v) if (!seen[v]) dfs(G, v);

        // back dfs
        int k = 0;
        groups.clear();
        seen.assign((int)G.size(), false);
        cmp.assign((int)G.size(), -1);
        for (int i = (int)G.size()-1; i >= 0; --i) {
            if (!seen[vs[i]]) {
                rvs.clear();
                rdfs(G, vs[i], k++);
                groups.push_back(rvs);
            }
        }
        reconstruct(G);
    }
};



//------------------------------//
// Examples
//------------------------------//

void Yosupo_Strongly_Connected_Components() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    
    // 入力
    int N, M;
    cin >> N >> M;
    Graph<int> G(N);
    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        G.add_edge(u, v);
    }
    //cout << G << endl;
    
    // SCC (not build dag)
    SCC<int> scc(G);
    const auto &groups = scc.groups;
    
    // 出力
    cout << groups.size() << endl;
    for (const auto &list : groups) {
        cout << list.size();
        for (auto v : list) cout << " " << v;
        cout << endl;
    }
}

void AOJ_GRL_3_C() {
    // 入力
    int N, M;
    cin >> N >> M;
    Graph<int> G(N);
    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        G.add_edge(u, v);
    }
    
    // SCC
    SCC<int> scc(G);
    
    // クエリ
    int Q;
    cin >> Q;
    while (Q--) {
        int u, v;
        cin >> u >> v;
        if (scc.cmp[u] == scc.cmp[v]) cout << 1 << endl;
        else cout << 0 << endl;
    }
}


int main() {
    //Yosupo_Strongly_Connected_Components();
    AOJ_GRL_3_C();
}