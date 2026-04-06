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
//   ARC 039 D - 旅行会社高橋君（for 二重辺連結成分分解）
//     https://arc039.contest.atcoder.jp/tasks/arc039_d
//
//   天下一プログラマーコンテスト2015予選A D - ハシポン（for 二重辺連結成分分解した Bridge-Block 木の考察）
//     https://atcoder.jp/contests/tenka1-2015-quala/tasks/tenka1_2015_qualA_d
//
//   TTPC 2024 DIV1 A - Don't Detect Cycle（for 橋列挙）
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

    // getter, bridge-block tree to orignal graph（v: node-id of bridge-block tree)
    int get_size(int v) const {
        return groups[v].size();
    }
    vector<int> get_group(int v) const {
        return groups[v];
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
        G.add_bidirected_edge(a, b, i);
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
        original_graph.add_bidirected_edge(s, t, 1);
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

// 天下一プログラマーコンテスト2015予選A D - ハシポン
void Tenka_2015_A_D() {
    int N, M, a, b;
    cin >> N >> M;
    Graph<int> G(N);
    for (int i = 0; i < M; i++) {
        cin >> a >> b;
        G.add_bidirected_edge(a, b);
    }
    TwoEdgeConnectedComponentsDecomposition<int> tecc(G);
    auto tree = tecc.tree;

    if (tree.size() == 1) {
        cout << "IMPOSSIBLE" << endl;
        return;
    } else if (tree.size() == 2) {
        cout << 0 << endl;
        return;
    } else if (tree.size() == 3) {
        bool allone = true;
        for (auto g : tecc.groups) if (g.size() > 1) allone = false;
        if (allone) cout << "IMPOSSIBLE" << endl;
        else cout << 1 << endl;
        return;
    }

    int leaf_num = 0;
    bool exist_length_one_bridge = false;
    for (int v = 0; v < tree.size(); v++) {
        if (tree[v].size() == 1) {
            leaf_num++;
            int v2 = tree[v][0];
            if (tree[v2].size() > 2) exist_length_one_bridge = true;
        }
    }
    if (exist_length_one_bridge) {
        cout << leaf_num / 2 << endl;
    } else {
        cout << (leaf_num + 1) / 2 << endl;
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
            G.add_bidirected_edge(u, v, i);
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
                    G2.add_bidirected_edge(e.from, e.to, e.val);
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
    YosupoLibraryCheckerTwoEdgeConnectedComponents();
    //ARC_039_D();
    //Tenka_2015_A_D();
    //TTPC_2024_DIV1_A();
}