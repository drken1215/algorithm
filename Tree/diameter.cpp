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


// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
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
        
    void add_edge(int from, int to, T val = 1) {
        list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }
    
    void add_bidirected_edge(int from, int to, T val = 1) {
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

// find diameter of tree
template<class Weight = long long> pair<Weight, vector<Edge<Weight>>> calc_diameter
(const Graph<Weight> &G) {
    Weight length = 0;
    vector<Edge<Weight>> path;
    vector<Edge<Weight>> prev(G.size(), Edge<Weight>(-1, -1));

    auto dfs = [&](auto &&dfs, int v, int p, bool record = true) -> pair<int, Weight> {
        pair<int, Weight> res{v, 0};
        for (const auto &e : G[v]) {
            if (e.to == p) continue;
            auto tmp = dfs(dfs, e.to, v, record);
            tmp.second += e.val;
            if (tmp.second > res.second) {
                res = tmp;
                if (record) prev[e.to] = e;
            }
        }
        return res;
    };

    auto [leaf, distance] = dfs(dfs, 0, -1, false);
    prev.assign((int)G.size(), Edge<Weight>(-1, -1));
    auto [most_distant_v, distance2] = dfs(dfs, leaf, -1, true);
    int cur = most_distant_v;
    while (cur != -1) {
        const auto &e = prev[cur];
        if (e.from == -1) break;
        length += e.val, path.emplace_back(e);
        cur = e.from;
    }
    reverse(path.begin(), path.end());
    return {length, path};
}


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
    auto [len, path] = calc_diameter(G);
    cout << len << " " << (int)path.size()+1 << endl;
    for (auto e : path) cout << e.from << " ";
    cout << path.back().to << endl;
}

// AGC 033 C - Removing Coins
void AGC_033_C() {
    int N, a, b;
    cin >> N;
    Graph G(N);
    for (int i = 0; i < N-1; ++i) {
        cin >> a >> b; --a, --b;
        G.add_bidirected_edge(a, b);
    }
    auto [len, path] = calc_diameter(G);   
    if (len % 3 == 1) cout << "Second" << endl;
    else cout << "First" << endl;
}


int main() {
    //Yosupo_Tree_Diameter();
    AGC_033_C();
}