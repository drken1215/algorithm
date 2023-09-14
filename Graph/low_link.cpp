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


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


using Graph = vector<vector<int> >;
struct LowLink {
    // main results
    vector<int> aps; // articulation points
    vector<pair<int,int> > brs; // brideges
    
    // intermediate results
    vector<int> seen, ord, low;
    void dfs_lowlink(const Graph &G, int v, int p = -1) {
        static int time = 0;
        seen[v] = true;
        ord[v] = low[v] = time++;
        int num_of_child = 0;
        bool exist = false; // for articulation point
        for (auto ch : G[v]) {
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
    void solve(const Graph &G) {
        int N = (int)G.size();
        seen.assign(N, 0); ord.resize(N); low.resize(N);
        aps.clear(); brs.clear();
        for (int v = 0; v < N; ++v) if (!seen[v]) dfs_lowlink(G, v);
    }
};


int main() {
    int V, E; cin >> V >> E;
    Graph G(V);
    for (int i = 0; i < E; ++i) {
        int s, t; cin >> s >> t;
        G[s].push_back(t);
        G[t].push_back(s);
    }
    LowLink ll;
    ll.solve(G);
    
    // 橋
    for (auto &br : ll.brs) if (br.first > br.second) swap(br.first, br.second);
    sort(ll.brs.begin(), ll.brs.end());
    for (auto br : ll.brs) cout << br.first << " " << br.second << endl;
    
    // 関節点
    sort(ll.aps.begin(), ll.aps.end());
    //for (auto v : ll.aps) cout << v << endl;
}
