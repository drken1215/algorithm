//
// サイクル検出
//
// cf.
//   DFS (深さ優先探索) 超入門！ 〜 グラフ・アルゴリズムの世界への入口 〜【後編】
//     https://qiita.com/drken/items/a803d4fc4a727e02f7ba
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
    vector<vector<Edge<T>>> reversed_list;  // reverse edge
    
    Graph(int n = 0) : list(n), reversed_list(n) { }
    void init(int n = 0) {
        list.assign(n, vector<Edge<T>>());
        reversed_list.assign(n, vector<Edge<T>>());
    }
    inline const vector<Edge<T>> &operator [] (int i) { return list[i]; }
    inline const vector<Edge<T>> &redges(int i) { return reversed_list[i]; }
    inline const size_t size() const { return list.size(); }
        
    void add_edge(int from, int to, T val = -1) {
        list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }
    
    void add_biedge(int from, int to, T val = -1) {
        list[from].push_back(Edge(from, to, val));
        list[to].push_back(Edge(to, from, val));
    }

    friend ostream &operator << (ostream &s, const Graph &G) {
        s << endl;
        for (int i = 0; i < G.size(); ++i) {
            s << i << ": " << G.list[i] << endl;
        }
        return s;
    }
};

// cycle detection
template<class T> struct CycleDetection {
    // input
    Graph<T> G;
    
    // intermediate results
    vector<bool> seen, finished;
    vector<Edge<T>> history;
    
    // constructor
    CycleDetection() { }
    CycleDetection(const Graph<T> &graph) { init(graph); }
    void init(const Graph<T> &graph) { G = graph; }
    
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
                return e2.to;
            }

            history.push_back(e2);
            int pos = dfs(e2.to, e2, is_prohibit_reverse);
            if (pos != -1) return pos;
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
    
    // find cycle
    vector<Edge<T>> detect(bool is_prohibit_reverse = true) {
        seen.assign(G.size(), false);
        finished.assign(G.size(), false);
        history.clear();
        int pos = -1;
        for (int v = 0; v < (int)G.size() && pos == -1; ++v) {
            if (seen[v]) continue;
            pos = dfs(v, Edge<T>(), is_prohibit_reverse);
            if (pos != -1) return reconstruct(pos);
        }
        return vector<Edge<T>>();
    }
};


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

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
        
        G.add_biedge(u, v, i);
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

void AOJ2891() {
    // 頂点数 (サイクルを一つ含むグラフなので辺数は N で確定)
    int N;
    cin >> N;

    // グラフ入力受取
    Graph<int> G(N);
    for (int i = 0; i < N; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        G.add_biedge(a, b, i);
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


int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    
    YosupoCycleDetectionDirected();
    //YosupoCycleDetectionUndirected();
    //AOJ2891();
}

