//
// 木の直径を求める
//
// verified:
//   Yosupo Library Checker - Tree Diameter
//     https://judge.yosupo.jp/problem/tree_diameter
//
//   AGC 033 C - Removing Coins
//     https://atcoder.jp/contests/agc033/tasks/agc033_c
//


#include <bits/stdc++.h>
using namespace std;


// find diameter of graph (not weighted)
struct Diameter {
    vector<int> path, prev;

    Diameter() {}
    Diameter(const vector<vector<int>> &G) {
        solve(G);
    }
    pair<int, int> DiameterDFS(const vector<vector<int>> &G, int v, int p) {
        pair<int, int> res(v, 0);
        for (auto to : G[v]) {
            if (to == p) continue;
            pair<int, int> tmp = DiameterDFS(G, to, v);
            tmp.second++;
            if (tmp.second > res.second) res = tmp, prev[to] = v;
        }
        return res;
    }
    vector<int> solve(const vector<vector<int>> &G) {
        prev.assign((int)G.size(), -1);
        auto [leaf, distance] = DiameterDFS(G, 0, -1);
        prev.assign((int)G.size(), -1);
        auto [ev, distance2] = DiameterDFS(G, leaf, -1);
        path.clear();
        int cur = ev;
        while (cur != -1) path.push_back(cur), cur = prev[cur];
        return path;
    }
};

// find diameter of weighted graph
// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = -1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
template<class T = long long> struct Graph {
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    
    Graph(int n = 0) : list(n), reversed_list(n) { }
    void init(int n = 0) {
        list.assign(n, vector<Edge<T>>());
        reversed_list.assign(n, vector<Edge<T>>());
    }
    Graph &operator = (const Graph &g) {
        list = g.list, reversed_list = g.reversed_list;
        return *this;
    }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    const vector<Edge<T>> &get_rev_edges(int i) const { return reversed_list[i]; }
    const size_t size() const { return list.size(); }
    const void clear() { list.clear(); }
    const void resize(int n) { list.resize(n); }
        
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

template<class T = long long> struct WeightedDiameter {
    vector<int> path;
    vector<pair<int, T>> prev;

    WeightedDiameter() {}
    WeightedDiameter(const Graph<T> &G) {
        solve(G);
    }
    pair<int, T> DiameterDFS(const Graph<T> &G, int v, int p) {
        pair<int, T> res{v, 0};
        for (auto e : G[v]) {
            if (e.to == p) continue;
            pair<int, T> tmp = DiameterDFS(G, e.to, v);
            tmp.second += e.val;
            if (tmp.second > res.second) res = tmp, prev[e.to] = {v, e.val};
        }
        return res;
    }
    pair<T, vector<int>> solve(const Graph<T> &G) {
        T res = 0;
        prev.assign((int)G.size(), make_pair(-1, -1));
        auto [leaf, distance] = DiameterDFS(G, 0, -1);
        prev.assign((int)G.size(), make_pair(-1, -1));
        auto [ev, distance2] = DiameterDFS(G, leaf, -1);
        path.clear();
        int cur = ev;
        while (cur != -1) {
            if (prev[cur].first != -1) res += prev[cur].second;
            path.push_back(cur), cur = prev[cur].first;
        }
        return {res, path};
    }
};


//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Tree Diameter
void Yosupo_Tree_Diameter() {
    int N, a, b;
    long long c;
    cin >> N;
    Graph<long long> G(N);
    for (int i = 0; i < N-1; i++) {
        cin >> a >> b >> c;
        G.add_bidirected_edge(a, b, c);
    }
    WeightedDiameter<long long> di;
    auto [len, path] = di.solve(G);
    cout << len << " " << (int)path.size() << endl;
    for (auto v : path) cout << v << " ";
    cout << endl;
}

// AGC 033 C - Removing Coins
void AGC_033_C() {
    int N, a, b;
    cin >> N;
    vector<vector<int>> G(N);
    for (int i = 0; i < N-1; ++i) {
        cin >> a >> b; --a, --b;
        G[a].push_back(b), G[b].push_back(a);
    }
    Diameter di;
    auto path = di.solve(G);    
    if (path.size() % 3 == 2) cout << "Second" << endl;
    else cout << "First" << endl;
}


int main() {
    Yosupo_Tree_Diameter();
    //AGC_033_C();
}