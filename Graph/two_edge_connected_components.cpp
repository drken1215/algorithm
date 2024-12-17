//
// Low-Link を用いた橋列挙と、二重辺連結成分分解
//
// cf.
//   hos: グラフ探索アルゴリズムとその応用
//     http://hos.ac/slides/20110504_graph.pdf
//
// verified:
//   ARC 039 D - 旅行会社高橋君
//     https://arc039.contest.atcoder.jp/tasks/arc039_d
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
    const void clear() {
        list.clear(), reversed_list.clear();
    }
    const void resize(int n) {
        list.resize(n), reversed_list.resize(n);
    }
        
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

// bridge decomposition
template<class T> struct BridgeDecomposition {
    // bridge and articulation points
    vector<int> aps; // articulation points
    vector<pair<int,int>> brs; // brideges
    
    // decomposition
    vector<vector<int>> scc; // scc[i] := i'th component
    vector<int> cmp;          // cmp[v] := which component is v belong
    Graph<T> newG;               // the tree
    
    // intermediate results
    vector<int> seen, ord, low;

    // constructor
    BridgeDecomposition() { }
    BridgeDecomposition(const Graph<T> &graph) {
        solve(graph);
    }
    void init(const Graph<T> &graph) {
        solve(graph);
    }
    
    // calc low-link
    void dfs_lowlink(const Graph<T> &G, int v, int p = -1) {
        static int time = 0;
        seen[v] = true;
        ord[v] = low[v] = time++;
        int num_of_child = 0;
        bool exist = false; // for articulation point
        for (auto ch_edge : G[v]) {
            int ch = ch_edge.to;
            if (seen[ch]) {
                if (ch != p) low[v] = min(low[v], ord[ch]); // back edge
                continue;
            }
            dfs_lowlink(G, ch, v);
            low[v] = min(low[v], low[ch]); // forward edge of DFS-tree
            if (ord[v] < low[ch]) brs.emplace_back(v, ch);
            if (ord[v] <= low[ch]) exist = true;
            ++num_of_child;
        }
        if ((p == -1 && num_of_child > 1) || (p != -1 && exist)) aps.emplace_back(v);
    }
    
    // reconstruct tree
    vector<int> comp;
    set<pair<int,int>> newEdges;
    void MakeTreeDFS(const Graph<T> &G, int v, int curcmp) {
        seen[v] = true;
        cmp[v] = curcmp;
        for (auto e : G[v]) {
            bool sameComp = true;
            if (seen[e.to]) sameComp = false;
            if (binary_search(brs.begin(), brs.end(), make_pair(v, e.to))
                || binary_search(brs.begin(), brs.end(), make_pair(e.to, v)))
                sameComp = false;
            
            if (!sameComp) {
                int newcmp = cmp[e.to];
                if (newcmp == -1) continue;
                if (newcmp != curcmp) {
                    if (!newEdges.count(make_pair(curcmp, newcmp))
                        && !newEdges.count(make_pair(newcmp, curcmp))) {
                        newEdges.insert(make_pair(curcmp, newcmp));
                        newG.add_bidirected_edge(curcmp, newcmp, e.val);
                    }
                }
            }
            else MakeTreeDFS(G, e.to, curcmp);
        }
        comp.push_back(v);
    }
    
    // main
    void solve(const Graph<T> &G) {
        // calc low-link
        int N = (int)G.size();
        seen.assign(N, 0), ord.resize(N), low.resize(N);
        aps.clear(), brs.clear();
        for (int v = 0; v < N; ++v) if (!seen[v]) dfs_lowlink(G, v);
        sort(brs.begin(), brs.end());
        
        // reconstruct tree
        scc.clear();
        newG.clear(), newG.resize(N);
        newEdges.clear();
        seen.assign(N, 0), cmp.assign(N, -1);
        int tcmp = 0;
        for (int v = 0; v < N; ++v) {
            if (seen[v]) continue;
            comp.clear();
            MakeTreeDFS(G, v, tcmp++);
            scc.push_back(comp);
        }
        newG.resize(tcmp);
    }
};



//------------------------------//
// Examples
//------------------------------//

// ARC 039 D - 旅行会社高橋君
void ARC_039_D() {
    int V, E, s, t;
    cin >> V >> E;
    Graph<int> original_graph(V);
    for (int i = 0; i < E; ++i) {
        cin >> s >> t, s--, t--;
        original_graph.add_bidirected_edge(s, t, 1);
    }
    BridgeDecomposition bd(original_graph);

    // LCA
    auto G = bd.newG;
    int h = 1;
    while ((1<<h) < V) ++h;
    vector<vector<int>> parent(h, vector<int>(V, -1));
    vector<int> depth(V, -1);

    auto dfs = [&](auto dfs, int v, int p, int d) -> void {
        parent[0][v] = p;
        depth[v] = d;
        for (auto e : G[v]) if (e.to != p) dfs(dfs, e.to, v, d+1);
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
        a = bd.cmp[a], b = bd.cmp[b], c = bd.cmp[c];
        int ab = dist(a, b);
        int bc = dist(b, c);
        int ac = dist(a, c);
        if (ab + bc == ac) puts("OK");
        else puts("NG");
    }
}


int main() {
    ARC_039_D();
}