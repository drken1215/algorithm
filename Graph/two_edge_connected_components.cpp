//
// Low-Link を用いた橋列挙と、二重辺連結成分分解
//
// cf.
//   hos: グラフ探索アルゴリズムとその応用
//     http://hos.ac/slides/20110504_graph.pdf
//
// verified:
//   Yosupo Library Checker - Two-Edge-Connected Components
//     https://judge.yosupo.jp/problem/two_edge_connected_components
//
//   ARC 039 D - 旅行会社高橋君
//     https://arc039.contest.atcoder.jp/tasks/arc039_d
//
//   TTPC 2024 DIV1 A - Don't Detect Cycle
//     https://atcoder.jp/contests/ttpc2024_1/tasks/ttpc2024_1_a
//


/*
    アイディア: DFS をしたとき、DFS 後退辺は橋とはなりえない
 
    ・ord[v] := 頂点を訪れた順番
    ・low[v] := v から「DFS 木の根から葉へ進む」or「後退辺を葉から根へ進む」ことによって辿り着ける頂点の ord の最小値
 
    DFS で u -> ... -> v と来て、v から u への後退辺があると、このサイクルの low がすべて ord[u] (以下) になる感じ
    このことから、
 
        DFS-search で、辺 v - ch を v -> ch の順に探索したときに、
            辺 v-to が橋　⇔　ord[v] < low[ch]
 */


#include <bits/stdc++.h>
using namespace std;


// Edge Class
template<class T> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1), val(-1) { }
    Edge(int f, int t, T v = -1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to;
    }
};

// graph class
template<class T> struct Graph {
    vector<vector<Edge<T>>> list;
    
    Graph() { }
    Graph(int n) : list(n) { }
    Graph(const Graph<T> &G) : list(G.list) { }
    void init(int n) {
        list.assign(n, vector<Edge<T>>());
    }
    Graph &operator = (const Graph &g) {
        list = g.list;
        return *this;
    }

    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    const size_t size() const { return list.size(); }
    const void clear() {
        list.clear();
    }
    const void resize(int n) {
        list.resize(n);
    }

    void add_edge(int from, int to, T val = -1) {
        list[from].push_back(Edge(from, to, val));
        list[to].push_back(Edge(to, from, val));
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

// low-link
template<class T> struct LowLink {
    // results
    vector<int> ord, low;
    vector<int> aps;         // articulation points
    vector<Edge<T>> brs;     // brideges

    // constructor
    LowLink() { }
    LowLink(const Graph<T> &G) {
        solve(G);
    }
    void init(const Graph<T> &G) {
        solve(G);
    }

    // solver
    int dfs(const Graph<T> &G, int t, int v, int p) {
        ord[v] = low[v] = t++;
        int num_of_children = 0;
        bool exist_articulation = false, is_multiple_edge = false;
        for (const auto &e : G[v]) {
            if (ord[e.to] == -1) {
                num_of_children++;
                t = dfs(G, t, e.to, v);
                low[v] = min(low[v], low[e.to]);  // forward edge of DFS-tree
                exist_articulation |= (p != -1) && (low[e.to] >= ord[v]);
                if (ord[v] < low[e.to]) brs.push_back(e);
            } else if (e.to != p || is_multiple_edge) {
                low[v] = min(low[v], ord[e.to]);  // back edge
            } else {
                is_multiple_edge = true;
            }
        }
        if (exist_articulation || (p == -1 && num_of_children > 1)) {
            aps.emplace_back(v);
        }
        return t;
    }

    void solve(const Graph<T> &G) {
        ord.assign(G.size(), -1), low.assign(G.size(), -1);
        for (int v = 0, k = 0; v < (int)G.size(); v++) {
            if (ord[v] == -1) k = dfs(G, k, v, -1);
        }
    }
};

// Two-Edge-Connected Components decomposition
template<class T> struct TwoEdgeConnectedComponentsDecomposition {
    // results
    LowLink<T> ll;
    vector<int> cmp;
    vector<vector<int>> groups, tree;

    // constructor
    TwoEdgeConnectedComponentsDecomposition() { }
    TwoEdgeConnectedComponentsDecomposition(const Graph<T> &G) {
        solve(G);
    }
    void init(const Graph<T> &G) {
        solve(G);
    }

    // solver
    int dfs(const Graph<T> &G, int t, int v, int p) {
        if (p >= 0 && ll.ord[p] >= ll.low[v]) cmp[v] = cmp[p];
        else cmp[v] = t++;
        for (const auto &e : G[v]) {
            if (cmp[e.to] == -1) t = dfs(G, t, e.to, v);
        }
        return t;
    }

    void solve(const Graph<T> &G) {
        ll.init(G);
        cmp.assign(G.size(), -1);
        int t = 0;
        for (int v = 0; v < (int)G.size(); v++) {
            if (cmp[v] == -1) t = dfs(G, t, v, -1);
        }
        groups.resize(t);
        tree.resize(t);
        for (int v = 0; v < (int)G.size(); v++) {
            groups[cmp[v]].push_back(v);
        }
        for (const auto &e : ll.brs) {
            int u = cmp[e.from], v = cmp[e.to];
            tree[u].push_back(v);
            tree[v].push_back(u);
        }
    }
};



//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Two-Edge-Connected Components
void YosupoLibraryCheckerTwoEdgeConnectedComponents() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N, M, a, b;
    cin >> N >> M;
    Graph<int> G(N);
    for (int i = 0; i < M; i++) {
        cin >> a >> b;
        G.add_edge(a, b, i);
    }
    TwoEdgeConnectedComponentsDecomposition<int> tecc(G);

    cout << tecc.groups.size() << '\n';
    for (const auto &group : tecc.groups) {
        cout << group.size();
        for (auto v : group) cout << " " << v;
        cout << '\n';
    }
}

// ARC 039 D - 旅行会社高橋君
void ARC_039_D() {
    int V, E, s, t;
    cin >> V >> E;
    Graph<int> original_graph(V);
    for (int i = 0; i < E; ++i) {
        cin >> s >> t, s--, t--;
        original_graph.add_edge(s, t, 1);
    }
    TwoEdgeConnectedComponentsDecomposition<int> tecc(original_graph);

    // LCA
    auto tree = tecc.tree;
    int h = 1;
    while ((1<<h) < V) ++h;
    vector<vector<int>> parent(h, vector<int>(V, -1));
    vector<int> depth(V, -1);

    auto dfs = [&](auto dfs, int v, int p, int d) -> void {
        parent[0][v] = p;
        depth[v] = d;
        for (auto to : tree[v]) if (to != p) dfs(dfs, to, v, d+1);
    };

    dfs(dfs, 0, -1, 0);
    for (int i = 0; i+1 < (int)parent.size(); ++i)
        for (int v = 0; v < V; ++v)
            if (parent[i][v] != -1)
                parent[i+1][v] = parent[i][parent[i][v]];

    auto get = [&](int u, int v) -> int {
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
    };

    auto dist = [&](int u, int v) -> int {
        int lca = get(u, v);
        return abs(depth[u] - depth[lca]) + abs(depth[v] - depth[lca]);
    };
    
    // queries
    int Q, a, b, c;
    cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        cin >> a >> b >> c, a--, b--, c--;
        a = tecc.cmp[a], b = tecc.cmp[b], c = tecc.cmp[c];
        int ab = dist(a, b);
        int bc = dist(b, c);
        int ac = dist(a, c);
        if (ab + bc == ac) puts("OK");
        else puts("NG");
    }
}

// TTPC 2024 DIV1 A - Don't Detect Cycle
void TTPC_2024_DIV1_A() {
    int T;
    cin >> T;

    while (T--) {
        int N, M, u, v;
        cin >> N >> M;
        Graph<int> G(N);
        for (int i = 0; i < M; i++) {
            cin >> u >> v, u--, v--;
            G.add_edge(u, v, i);
        }

        vector<int> res, former, latter;
        while (true) {
            set<pair<int,int>> ers;

            // erase bridges
            LowLink<int> ll(G);
            for (auto e : ll.brs) {
                former.push_back(e.val);
                ers.insert(minmax(e.from, e.to));
            }

            // erase edges whose endpoints with 2 degree
            vector<int> deg(N, 0);
            for (int v = 0; v < N; v++) {
                for (auto e : G[v]) {
                    if (ers.count(minmax(e.from, e.to))) continue;
                    deg[v]++;
                }
            }
            for (int v = 0; v < N; v++) {
                for (auto e : G[v]) {
                    if (e.from > e.to) continue;
                    if (ers.count(minmax(e.from, e.to))) continue;
                    if (deg[e.from] == 2 && deg[e.to] == 2) {
                        latter.push_back(e.val);
                        ers.insert(minmax(e.from, e.to));
                    }
                }
            }

            if (ers.empty()) break;

            // build new graph
            Graph<int> G2(N);
            for (int v = 0; v < N; v++) {
                for (auto e : G[v]) {
                    if (e.from > e.to) continue;
                    if (ers.count(minmax(e.from, e.to))) continue;
                    G2.add_edge(e.from, e.to, e.val);
                }
            }
            G = G2;
        }

        reverse(latter.begin(), latter.end());
        for (auto v : former) res.push_back(v);
        for (auto v : latter) res.push_back(v);
        if (res.size() < M) {
            cout << -1 << endl;
        } else {
            for (auto id : res) cout << id+1 << " ";
            cout << endl;
        }
    }
}


int main() {
    //YosupoLibraryCheckerTwoEdgeConnectedComponents();
    //ARC_039_D();
    TTPC_2024_DIV1_A();
}