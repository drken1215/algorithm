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


#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
using namespace std;


#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, vector<vector<T> > P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }




using Edge = pair<int,long long>;
using Graph = vector<vector<Edge> >;

struct BridgeDecomposition {
    // bridge and articulation points
    vector<int> aps; // articulation points
    vector<pair<int,int> > brs; // brideges
    
    // decomposition
    vector<vector<int> > scc; // scc[i] := i'th component
    vector<int> cmp;          // cmp[v] := which component is v belong
    Graph newG;               // the tree
    
    // intermediate results
    vector<int> seen, ord, low;
    
    // calc low-link
    void dfs_lowlink(const Graph &G, int v, int p = -1) {
        static int time = 0;
        seen[v] = true;
        ord[v] = low[v] = time++;
        int num_of_child = 0;
        bool exist = false; // for articulation point
        for (auto ch_edge : G[v]) {
            int ch = ch_edge.first;
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
    set<pair<int,int> > newEdges;
    void MakeTreeDFS(const Graph &G, int v, int curcmp) {
        seen[v] = true;
        cmp[v] = curcmp;
        for (auto e : G[v]) {
            int to = e.first;
            auto weight = e.second;
            
            bool sameComp = true;
            if (seen[to]) sameComp = false;
            if (binary_search(brs.begin(), brs.end(), make_pair(v, to))
                || binary_search(brs.begin(), brs.end(), make_pair(to, v)))
                sameComp = false;
            
            if (!sameComp) {
                int newcmp = cmp[to];
                if (newcmp == -1) continue;
                if (newcmp != curcmp) {
                    if (!newEdges.count(make_pair(curcmp, newcmp))
                        && !newEdges.count(make_pair(newcmp, curcmp))) {
                        newEdges.insert(make_pair(curcmp, newcmp));
                        newG[curcmp].push_back(Edge(newcmp, weight));
                        newG[newcmp].push_back(Edge(curcmp, weight));
                    }
                }
            }
            else MakeTreeDFS(G, to, curcmp);
        }
        comp.push_back(v);
    }
    
    // main
    void solve(const Graph &G) {
        // calc low-link
        int N = (int)G.size();
        seen.assign(N, 0); ord.resize(N); low.resize(N);
        aps.clear(); brs.clear();
        for (int v = 0; v < N; ++v) if (!seen[v]) dfs_lowlink(G, v);
        sort(brs.begin(), brs.end());
        
        // reconstruct tree
        scc.clear();
        newG.clear(); newG.resize(N);
        newEdges.clear();
        seen.assign(N, 0); cmp.assign(N, -1);
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


// LCA
struct LCA {
    vector<vector<int> > parent; // parent[d][v] := 2^d-th parent of v
    vector<int> depth;
    LCA() { }
    LCA(const Graph &G, int r = 0) { init(G, r); }
    void init(const Graph &G, int r = 0) {
        int V = (int)G.size();
        int h = 1;
        while ((1<<h) < V) ++h;
        parent.assign(h, vector<int>(V, -1));
        depth.assign(V, -1);
        dfs(G, r, -1, 0);
        for (int i = 0; i+1 < (int)parent.size(); ++i)
            for (int v = 0; v < V; ++v)
                if (parent[i][v] != -1)
                    parent[i+1][v] = parent[i][parent[i][v]];
    }
    void dfs(const Graph &G, int v, int p, int d) {
        parent[0][v] = p;
        depth[v] = d;
        for (auto e : G[v]) if (e.first != p) dfs(G, e.first, v, d+1);
    }
    int get(int u, int v) {
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
    }
    int dist(int u, int v) {
        int lca = get(u, v);
        return abs(depth[u] - depth[lca]) + abs(depth[v] - depth[lca]);
    }
};


int main() {
    int V, E; cin >> V >> E;
    Graph G(V);
    for (int i = 0; i < E; ++i) {
        int s, t; cin >> s >> t; --s, --t;
        G[s].push_back(Edge(t, 1));
        G[t].push_back(Edge(s, 1));
    }
    
    BridgeDecomposition bd;
    bd.solve(G);
    
    /*
    COUT(bd.brs);
    COUT(bd.newG);
    COUT(bd.scc);
    COUT(bd.cmp);
    */
    
    LCA lca(bd.newG);
    int Q; cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        int a, b, c; cin >> a >> b >> c; --a, --b, --c;
        a = bd.cmp[a]; b = bd.cmp[b]; c = bd.cmp[c];
        int ab = lca.dist(a, b);
        int bc = lca.dist(b, c);
        int ac = lca.dist(a, c);
        if (ab + bc == ac) puts("OK");
        else puts("NG");
    }
}
