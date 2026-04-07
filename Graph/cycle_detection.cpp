//
// サイクル検出
//
// cf.
//   グラフのサイクル検出 (閉路検出) by DFS
//     https://drken1215.hatenablog.com/entry/2023/05/20/200517
//
// Verified:
//   Yosupo Library Checker - Cycle Detection (Directed)
//     https://judge.yosupo.jp/problem/cycle_detection
//
//   Yosupo Library Checker - Cycle Detection (Undirected)
//     https://judge.yosupo.jp/problem/cycle_detection_undirected
//
//   AOJ 2891 - な◯りカット
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2891
//
//   ABC 256 E - Takahashi's Anguish
//     https://atcoder.jp/contests/abc256/tasks/abc256_e
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
    int V, E;
    bool record_reversed_edges = false, record_edge_index = false;
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    vector<unordered_map<int, int>> id;  // id[v][w] := the index of node w in G[v]

    // constructors
    Graph(int n = 0, int m = 0, bool rre = false, bool rei = false) { 
        init(n, m, rre, rei);
    }
    void init(int n = 0, int m = 0, bool rre = false, bool rei = false) {
        V = n, E = m, record_reversed_edges = rre, record_edge_index = rei;
        list.assign(n, vector<Edge<T>>());
        if (record_reversed_edges) reversed_list.assign(n, vector<Edge<T>>());
        if (record_edge_index) id.assign(n, unordered_map<int, int>());
    }
    Graph(const Graph&) = default;
    Graph& operator = (const Graph&) = default;

    // getters
    vector<Edge<T>> &operator [] (int i) { return list[i]; }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    constexpr size_t size() const { return list.size(); }
    constexpr void clear() { V = 0; list.clear(); }
    constexpr void resize(int n) { V = n; list.resize(n); }
    const vector<Edge<T>> &get_rev_edges(int i) const { 
        assert(record_reversed_edges);
        return reversed_list[i];
    }
    Edge<T> &get_edge(int u, int v) {
        assert(record_edge_index);
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u][v] >= 0 && id[u][v] < list[u].size());
        return list[u][id[u][v]];
    }
    const Edge<T> &get_edge(int u, int v) const {
        assert(record_edge_index);
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u].at(v) >= 0 && id[u].at(v) < list[u].size());
        return list[u][id[u].at(v)];
    }

    // add edge
    void add_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size(); 
        list[from].push_back(Edge(from, to, val));
        if (record_reversed_edges) reversed_list[to].push_back(Edge(to, from, val));
    }
    void add_bidirected_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size(); 
        list[from].push_back(Edge(from, to, val));
        if (record_reversed_edges) reversed_list[from].push_back(Edge(from, to, val));
        if (from != to) {
            if (record_edge_index) id[to][from] = (int)list[to].size(); 
            list[to].push_back(Edge(to, from, val));
            if (record_reversed_edges) reversed_list[to].push_back(Edge(to, from, val));
        }
    }

    // input / output
    friend istream& operator >> (istream &is, Graph &G) {
        for (int i = 0; i < G.E; i++) {
            int u, v;
            is >> u >> v, u--, v--;
            G.add_bidirected_edge(u, v);
        }
        return is;
    }
    friend ostream &operator << (ostream &os, const Graph &G) {
        os << endl;
        for (int i = 0; i < (int)G.size(); ++i) {
            os << i << " -> ";
            for (int j = 0; j < (int)G[i].size(); j++) {
                if (j) os << ", ";
                os << G[i][j].to << "(" << G[i][j].val << ")";
            }
            os << endl;
        }
        return os;
    }
};

// cycle detection
template<class T = long long> struct CycleDetection {
    // input
    Graph<T> G;
    
    // intermediate results
    vector<bool> seen, finished;
    vector<Edge<T>> history;
    
    // constructor
    CycleDetection() { }
    CycleDetection(const Graph<T> &graph) { init(graph); }
    void init(const Graph<T> &graph) {
        G = graph;
        seen.assign(G.size(), false);
        finished.assign(G.size(), false);
    }
    
    // dfs
    // return the vertex where cycle is detected
    int dfs(int v, const Edge<T> &e, bool is_prohibit_reverse = true) {
        seen[v] = true;
        for (const Edge<T> &e2 : G[v]) {
            if (is_prohibit_reverse && e2.to == e.from) continue;
            if (finished[e2.to]) continue;

            // detect cycle
            if (seen[e2.to] && !finished[e2.to]) {
                history.push_back(e2);
                finished[v] = true;
                return e2.to;
            }

            history.push_back(e2);
            int pos = dfs(e2.to, e2, is_prohibit_reverse);
            if (pos != -1) {
                finished[v] = true;
                return pos;
            }
            history.pop_back();
        }
        finished[v] = true;
        return -1;
    }
    
    // reconstruct
    vector<Edge<T>> reconstruct(int pos) {
        vector<Edge<T>> cycle;
        while (!history.empty()) {
            const Edge<T> &e = history.back();
            cycle.push_back(e);
            history.pop_back();
            if (e.from == pos) break;
        }
        reverse(cycle.begin(), cycle.end());
        return cycle;
    }
    
    // find cycle, v is the start vertex
    vector<Edge<T>> detect_from_v(int v, bool is_prohibit_reverse = true) {
        history.clear();
        int pos = dfs(v, Edge<T>(), is_prohibit_reverse);
        if (pos != -1) return reconstruct(pos);
        else return vector<Edge<T>>();
    }
    
    // find cycle
    vector<Edge<T>> detect(bool is_prohibit_reverse = true) {
        int pos = -1;
        for (int v = 0; v < (int)G.size() && pos == -1; ++v) {
            if (seen[v]) continue;
            history.clear();
            pos = dfs(v, Edge<T>(), is_prohibit_reverse);
            if (pos != -1) return reconstruct(pos);
        }
        return vector<Edge<T>>();
    }
};



//------------------------------//
// Examples
//------------------------------//

void YosupoCycleDetectionDirected() {
    // 有向グラフの受け取り
    int N, M;
    cin >> N >> M;
    Graph<int> G(N);
    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        G.add_edge(u, v, i);
    }
    
    
    // cycle detection
    CycleDetection<int> cd(G);
    
    // false: accept (u, v), (v, u)
    const vector<Edge<int>> &res = cd.detect(false);
    
    // 出力
    if (res.empty()) cout << -1 << endl;
    else {
        cout << res.size() << endl;
        for (const Edge<int> &e : res) {
            cout << e.val << endl;
        }
    }
}

void YosupoCycleDetectionUndirected() {
    // 有向グラフの受け取り, 多重辺と自己ループは予め除去しておく
    int N, M;
    cin >> N >> M;
    Graph<int> G(N);
    map<pair<int,int>, int> edges;
    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        
        // self-loop
        if (u == v) {
            cout << 1 << endl << u << endl << i << endl;
            return;
        }
        
        // multi-edge
        if (u > v) swap(u, v);
        if (edges.count({u, v})) {
            cout << 2 << endl;
            cout << u << " " << v << endl;
            cout << edges[{u, v}] << " " << i << endl;
            return;
        }
        
        G.add_bidirected_edge(u, v, i);
        edges[{u, v}] = i;
    }
    
    // cycle detection
    CycleDetection<int> cd(G);
    const vector<Edge<int>> &res = cd.detect(true);
    
    // 出力
    if (res.empty()) cout << -1 << endl;
    else {
        cout << res.size() << endl;
        for (const Edge<int> &e : res) cout << e.from << " ";
        cout << endl;
        for (const Edge<int> &e : res) cout << e.val << " ";
        cout << endl;
    }
}

void AOJ_2891() {
    // 頂点数 (サイクルを一つ含むグラフなので辺数は N で確定)
    int N;
    cin >> N;

    // グラフ入力受取
    Graph<int> G(N);
    for (int i = 0; i < N; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        G.add_bidirected_edge(a, b, i);
    }

    // 探索
    CycleDetection<int> cd(G);
    const vector<Edge<int>> &cycle = cd.detect();
    set<int> namori;
    for (const Edge<int> &e : cycle) namori.insert(e.to);

    // クエリに答える
    int Q;
    cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        if (namori.count(a) && namori.count(b)) cout << 2 << endl;
        else cout << 1 << endl;
    }
}

void ABC_256_E() {
    int N;
    cin >> N;
    Graph<long long> G(N);
    vector<long long> X(N), C(N);
    for (int i = 0; i < N; ++i) cin >> X[i];
    for (int i = 0; i < N; ++i) cin >> C[i];
    for (int i = 0; i < N; ++i) G.add_edge(i, X[i]-1, C[i]);

    long long res = 0;
    CycleDetection<long long> cd(G);
    for (int v = 0; v < N; ++v) {
        if (cd.seen[v]) continue;
        const auto &cycle = cd.detect_from_v(v, false);
        if (cycle.empty()) continue;
        long long minv = 1LL<<60;
        for (const auto &e : cycle) minv = min(minv, e.val);
        res += minv;
    }
    cout << res << endl;
}


int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    
    YosupoCycleDetectionDirected();
    //YosupoCycleDetectionUndirected();
    //AOJ_2891();
    //ABC_256_E();
}

