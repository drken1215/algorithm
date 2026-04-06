//
// Low-Link を用いた、橋・関節点列挙, O(V + E)
//
// cf.
//   hos: グラフ探索アルゴリズムとその応用
//     http://hos.ac/slides/20110504_graph.pdf
//
// verified:
//   AOJ Course GRL_3_B Connected Components - Bridges
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_3_B&lang=jp
//
//   AOJ Course GRL_3_A Connected Components - Articulation Points
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_3_A&lang=jp
//
//


/*
    アイディア: DFS をしたとき、DFS 後退辺は橋とはなりえない
 
    ・ord[v] := 頂点を訪れた順番
    ・low[v] := v から「DFS 木の根から葉へ進む」or「後退辺を葉から根へ進む」ことによって辿り着ける頂点の ord の最小値
 
    DFS で u -> ... -> v と来て、v から u への後退辺があると、このサイクルの low がすべて ord[u] (以下) になる感じ
    このtことから、
 
        DFS-search で、辺 v - ch を v -> ch の順に探索したときに、
            辺 v-to が橋　⇔　ord[v] < low[ch]
 
        DFS-search で、
            頂点 v が関節点　⇔
                ・v が根のとき、deg[v] > 1
                .それ以外のとき、u のある子供 ch が存在して、ord[v] <= low[ch]
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



//------------------------------//
// Examples
//------------------------------//

// AOJ Course GRL_3_B Connected Components - Bridges
void AOJ_GRL_3_B() {
    int V, E;
    cin >> V >> E;
    Graph G(V);
    for (int i = 0; i < E; ++i) {
        int s, t;
        cin >> s >> t;
        G.add_bidirected_edge(s, t);
    }
    LowLink link(G);
    
    // 橋
    vector<pair<int,int>> res;
    for (auto &e : link.brs) res.emplace_back(min(e.from, e.to), max(e.from, e.to));
    sort(res.begin(), res.end());
    for (auto [u, v] : res) cout << u << " " << v << endl;
}

// AOJ Course GRL_3_A Connected Components - Articulation Points
void AOJ_GRL_3_A() {
    int V, E;
    cin >> V >> E;
    Graph G(V);
    for (int i = 0; i < E; ++i) {
        int s, t;
        cin >> s >> t;
        G.add_bidirected_edge(s, t);
    }
    LowLink link(G);
    
    // 関節点
    sort(link.aps.begin(), link.aps.end());
    for (auto v : link.aps) cout << v << endl;
}


int main() {
    AOJ_GRL_3_B();
    //AOJ_GRL_3_A();
}
