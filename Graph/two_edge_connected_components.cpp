//
// Low-Link を用いた橋列挙と、二重辺連結成分分解・二重頂点連結成分分解
//
// cf.
//   hos: グラフ探索アルゴリズムとその応用
//     http://hos.ac/slides/20110504_graph.pdf
//
// verified:
//   Yosupo Library Checker - Two-Edge-Connected Components
//     https://judge.yosupo.jp/problem/two_edge_connected_components
//
//   Yosupo Library Checker - Biconnected Components
//     https://judge.yosupo.jp/problem/biconnected_components
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
//   ARC 045 D - みんな仲良し高橋君（for 二重頂点連結成分分解した Block-Cut 木上の DP）
//     https://atcoder.jp/contests/arc045/tasks/arc045_d
//
//   AOJ 3022 Problem J: Cluster Network（for 二重頂点連結成分分解した Block-Cut 木上の DP）
//     https://onlinejudge.u-aizu.ac.jp/problems/3022
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

// BiConnected Components decomposition
// block-cut tree (aps: 0, 1, ..., A-1, components: A, A+1, ..., A+C-1)
// (A: size of aps, C: num of components)
template<class T> struct BiConnectedComponentsDecomposition {
    // result
    LowLink<T> ll;
    vector<int> id_ap;   // index of the articulation point (size: V)
    vector<int> id_cc;   // index of the connected component (size: V)
    vector<vector<int>> groups;   // biconnected components (size: C)
    vector<vector<int>> tree;     // block-cut tree (size: A + C)

    // intermediate results
    vector<int> seen, finished;
    vector<vector<pair<int, int>>> grouped_edges;
    vector<pair<int, int>> tmp_edges;

    // constructor
    BiConnectedComponentsDecomposition() { }
    BiConnectedComponentsDecomposition(const Graph<T> &G) {
        solve(G);
    }
    void init(const Graph<T> &G) {
        solve(G);
    }

    // getter, original graph to block-cut tree (v: node of orignal graph)
    int is_ap_original_graph(int v) const {
        return (id_ap[v] != -1);
    }
    int get_id(int v) const {
        return (id_ap[v] == -1 ? id_cc[v] : id_ap[v]);
    }

    // getter, block-cut tree to orignal graph（v: node-id of block-cut tree)
    int is_ap(int v) const {
        return (v < ll.aps.size());
    }
    int get_ap(int v) const {
        if (v >= (int)ll.aps.size()) return -1;  // not ap
        else return ll.aps[v];
    }
    int get_size(int v) const {  // including aps
        if (v < (int)ll.aps.size()) return 1;  // ap
        else return groups[v - ll.aps.size()].size();
    }
    vector<int> get_group(int v) const {
        if (v < (int)ll.aps.size()) return vector<int>({ll.aps[v]});  // ap
        else return groups[v - ll.aps.size()];
    }

    // solver
    void dfs(const Graph<T> &G, int v, int p) {
        seen[v] = true;
        if (G[v].empty()) {
            groups.emplace_back(vector<int>({v}));
        }
        for (const auto &e : G[v]) {
            if (e.to == p) continue;
            if (!seen[e.to] || ll.ord[e.to] < ll.ord[v]) {
                tmp_edges.emplace_back(minmax(v, e.to));
            }
            if (!seen[e.to]) {
                dfs(G, e.to, v);
                if (ll.low[e.to] >= ll.ord[v]) {
                    groups.emplace_back(vector<int>({v}));
                    grouped_edges.emplace_back();
                    int ap = v;
                    while (!tmp_edges.empty()) {
                        const auto &e2 = tmp_edges.back();
                        if (!finished[e2.first] && e2.first != ap) {
                            groups.back().emplace_back(e2.first);
                            finished[e2.first] = true;
                        }
                        if (!finished[e2.second] && e2.second != ap) {
                            groups.back().emplace_back(e2.second);
                            finished[e2.second] = true;
                        }
                        grouped_edges.back().emplace_back(e2);
                        tmp_edges.pop_back();
                        if (e2.first == min(v, e.to) && e2.second == max(v, e.to)) break;
                    }
                }
            }
        }
    }

    void solve(const Graph<T> &G) {
        ll.init(G);
        seen.assign(G.size(), false), finished.assign(G.size(), false);
        for (int v = 0; v < (int)G.size(); v++) {
            if (!seen[v]) dfs(G, v, -1);
        }
        id_ap.assign(G.size(), -1), id_cc.assign(G.size(), -1);
        for (int i = 0; i < (int)ll.aps.size(); i++) {
            id_ap[ll.aps[i]] = i;
        }
        tree.assign(ll.aps.size() + grouped_edges.size(), vector<int>());
        vector<int> last(G.size(), -1);
        for (int i = 0; i < (int)grouped_edges.size(); i++) {
            vector<int> st;
            for (auto [u, v] : grouped_edges[i]) {
                st.push_back(u), st.push_back(v);
            }
            for (auto v : st) {
                if (id_ap[v] == -1) {
                    id_cc[v] = i + ll.aps.size();
                } else if (last[v] != i) {
                    tree[i + ll.aps.size()].push_back(id_ap[v]);
                    tree[id_ap[v]].push_back(i + ll.aps.size());
                    last[v] = i;
                }
            }
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

// Yosupo Library Checker - BiConnected Components
void YosupoLibraryCheckerBiConnectedComponents() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N, M, a, b;
    cin >> N >> M;
    Graph<int> G(N);
    for (int i = 0; i < M; i++) {
        cin >> a >> b;
        G.add_edge(a, b, i);
    }
    BiConnectedComponentsDecomposition<int> bcc(G);

    cout << bcc.groups.size() << '\n';
    for (const auto &group : bcc.groups) {
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

// 天下一プログラマーコンテスト2015予選A D - ハシポン
void Tenka_2015_A_D() {
    int N, M, a, b;
    cin >> N >> M;
    Graph<int> G(N);
    for (int i = 0; i < M; i++) {
        cin >> a >> b;
        G.add_edge(a, b);
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


// Union-Find
struct UnionFind {
    // core member
    vector<int> par, nex;

    // constructor
    UnionFind() { }
    UnionFind(int N) : par(N, -1), nex(N) {
        init(N);
    }
    void init(int N) {
        par.assign(N, -1);
        nex.resize(N);
        for (int i = 0; i < N; ++i) nex[i] = i;
    }
    
    // core methods
    int root(int x) {
        if (par[x] < 0) return x;
        else return par[x] = root(par[x]);
    }
    
    bool same(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y) {
        x = root(x), y = root(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y); // merge technique
        par[x] += par[y];
        par[y] = x;
        swap(nex[x], nex[y]);
        return true;
    }
    
    int size(int x) {
        return -par[root(x)];
    }
    
    // get group
    vector<int> group(int x) {
        vector<int> res({x});
        while (nex[res.back()] != x) res.push_back(nex[res.back()]);
        return res;
    }
    vector<vector<int>> groups() {
        vector<vector<int>> member(par.size());
        for (int v = 0; v < (int)par.size(); ++v) {
            member[root(v)].push_back(v);
        }
        vector<vector<int>> res;
        for (int v = 0; v < (int)par.size(); ++v) {
            if (!member[v].empty()) res.push_back(member[v]);
        }
        return res;
    }
    
    // debug
    friend ostream& operator << (ostream &s, UnionFind uf) {
        const vector<vector<int>> &gs = uf.groups();
        for (const vector<int> &g : gs) {
            s << "group: ";
            for (int v : g) s << v << " ";
            s << endl;
        }
        return s;
    }
};

void ARC_045_D() {
    // 入力
    int N;
    cin >> N;
    vector<int> X(N*2+1), Y(N*2+1);
    for (int i = 0; i < N*2+1; i++) cin >> X[i] >> Y[i];

    // 横方向
    using pint = pair<int,int>;
    vector<vector<pint>> xy(N*2+2), yx(N*2+2);
    for (int i = 0; i < N*2+1; i++) {
        xy[X[i]].emplace_back(Y[i], i);
        yx[Y[i]].emplace_back(X[i], i);
    }
    for (int v = 0; v <= N*2+1; v++) {
        sort(xy[v].begin(), xy[v].end());
        sort(yx[v].begin(), yx[v].end());
    }

    // グラフを作る（横・縦に「隣接」「1個飛ばし」のみ辺を張る）
    Graph<int> G(N*2+1);
    for (int v = 0; v <= N*2+1; v++) {
        for (int i = 0; i < xy[v].size(); i++) {
            if (i+1 < xy[v].size()) G.add_edge(xy[v][i].second, xy[v][i+1].second);
            if (i+2 < xy[v].size()) G.add_edge(xy[v][i].second, xy[v][i+2].second);
        }
        for (int i = 0; i < yx[v].size(); i++) {
            if (i+1 < yx[v].size()) G.add_edge(yx[v][i].second, yx[v][i+1].second);
            if (i+2 < yx[v].size()) G.add_edge(yx[v][i].second, yx[v][i+2].second);
        }
    }

    // Block-Cut 木の構築
    BiConnectedComponentsDecomposition<int> bcc(G);

    // グラフ G の各連結成分のサイズを求めて、間接点以外の点の答えを求める
    UnionFind uf(N*2+1);
    for (int v = 0; v < N*2+1; v++) for (auto e : G[v]) uf.merge(e.from, e.to);
    int odd_num = 0;
    for (int v = 0; v < N*2+1; v++) if (uf.root(v) == v && uf.size(v) % 2 == 1) odd_num++;
    if (odd_num > 1) {
        // 奇数サイズの連結成分が複数個ある場合はすべて NG（以降、奇数サイズの連結成分は 1 個とする)
        for (int v = 0; v < N*2+1; v++) cout << "NG" << endl;
        return;
    }
    vector<bool> res(N*2+1, false);
    for (int v = 0; v < N*2+1; v++) {
        // 奇数サイズかつ間接点以外は true (間接点の結果はあとで上書きする)
        if (uf.size(v) % 2 == 1) {
            res[v] = true;
        }
    }

    // Block-Cut 木上の探索により、間接点の答えを求める
    auto tree = bcc.tree;
    vector<bool> seen(tree.size(), false);
    auto dfs = [&](auto dfs, int v, int p) -> int {
        seen[v] = true;

        // 奇数サイズの連結成分を二重頂点連結成分分解してできる Block-Cut 木を
        // 根付き木として探索したとき、各間接点について、
        // ある子頂点が存在して、それを根とする部分木のサイズが奇数ならば、"NG"
        int siz = bcc.get_size(v), odd_num = 0;
        for (auto ch : tree[v]) {
            if (ch == p) continue;
            int tmp = dfs(dfs, ch, v) - 1;
            if (tmp % 2 == 1) odd_num++;
            siz += tmp;
        }
        if (bcc.is_ap(v)) {
            int apv = bcc.get_ap(v);  // Block-cut 木の頂点 v に対応するもとのグラフの頂点
            if (uf.size(apv) % 2 == 1 && odd_num > 0) {
                res[apv] = false;
            }
        }
        return siz;
    };
    for (int v = 0; v < tree.size(); v++) {
        if (seen[v]) continue;
        dfs(dfs, v, -1);
    } 

    // 出力
    for (int v = 0; v < res.size(); v++) {
        if (res[v]) cout << "OK" << endl;
        else cout << "NG" << endl;
    }
}

// AOJ 3022 - Cluster Network
void AOJ_3022() {
    long long N, M, u, v, all = 0;
    cin >> N >> M;
    vector<long long> w(N), res(N, 0);
    for (int i = 0; i < N; i++) cin >> w[i], all += w[i];
    Graph<int> G(N);
    for (int i = 0; i < M; i++) {
        cin >> u >> v, u--, v--;
        G.add_edge(u, v), G.add_edge(v, u);
    }

    // Block-Cut 木上を形成
    BiConnectedComponentsDecomposition<int> bcc(G);

    // 関節点以外について求める
    for (int ov = 0; ov < N; ov++) {
        if (!bcc.is_ap_original_graph(ov)) res[ov] = all - w[ov];
    }

    // 関節点について求める：Block-Cut 木上の DP
    auto tree = bcc.tree;
    vector<long long> sum(tree.size(), 0);
    for (int v = 0; v < tree.size(); v++) {
        const auto &group = bcc.get_group(v);
        for (auto ov : group) sum[v] += w[ov];
    }
    auto rec = [&](auto rec, int v, int p) -> long long {
        const auto &group = bcc.get_group(v);
        long long ma = 0, all_weight = sum[v];
        for (auto ch : tree[v]) {
            if (ch == p) continue;
            long long sub = rec(rec, ch, v);
            if (bcc.is_ap(v)) sub -= sum[v];
            else sub -= sum[ch];
            ma = max(ma, sub), all_weight += sub;
        }

        // 関節点について処理する
        if (bcc.is_ap(v)) {
            int ov = bcc.get_ap(v);
            long long rem = all - all_weight;
            ma = max(ma, rem);
            res[ov] = ma;
        }
        return all_weight;
    };
    rec(rec, 0, -1);

    for (auto val : res) cout << val << endl;
}

int main() {
    //YosupoLibraryCheckerTwoEdgeConnectedComponents();
    //YosupoLibraryCheckerBiConnectedComponents();
    //ARC_039_D();
    Tenka_2015_A_D();
    //TTPC_2024_DIV1_A();
    //ARC_045_D();
    //AOJ_3022();
}