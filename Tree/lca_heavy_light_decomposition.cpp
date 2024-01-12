//
// LCA (by HL 分解)
//
// verified:
//   AOJ Course GRL_5_C Tree - Lowest Common Ancestor
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_5_C&lang=jp
//


#include <bits/stdc++.h>
using namespace std;


// HL-Decomposition
// vid: id of v after HL-Decomposition
// inv: inv[vid[v]] = v
// par: id of parent
// depth
// subsize: size of subtree
// head: head-id in the heavy-path
// prev, next: prev-id, next-id in the heavy-path
// type: the id of tree for forest
// vend: the last-id of node in v-subtree

using Graph = vector<vector<int>>;
struct HLDecomposition {
    int n;
    Graph G;
    vector<int> vid, inv, par, depth, subsize, head, prev, next, type;
    
    // construct
    HLDecomposition() { }
    HLDecomposition(const Graph &G_) :
        n((int)G_.size()), G(G_),
        vid(n, -1), inv(n), par(n), depth(n), subsize(n, 1),
        head(n), prev(n, -1), next(n, -1), type(n) { }
    void build(vector<int> roots = {0}) {
        int curtype = 0, pos = 0;
        for (auto r : roots) decide_heavy_edge(r), reconstruct(r, curtype++, pos);
    }
    void decide_heavy_edge(int r) {
        stack<pair<int,int> > st;
        par[r] = -1, depth[r] = 0;
        st.emplace(r, 0);
        while (!st.empty()) {
            int v = st.top().first;
            int &i = st.top().second;
            if (i < (int)G[v].size()) {
                int e = G[v][i++];
                if (e == par[v]) continue;
                par[e] = v, depth[e] = depth[v] + 1;
                st.emplace(e, 0);
            }
            else {
                st.pop();
                int maxsize = 0;
                for (auto e : G[v]) {
                    if (e == par[v]) continue;
                    subsize[v] += subsize[e];
                    if (maxsize < subsize[e]) maxsize = subsize[e], prev[e] = v, next[v] = e;
                }
            }
        }
    }
    void reconstruct(int r, int curtype, int &pos) {
        stack<int> st({r});
        while (!st.empty()) {
            int start = st.top(); st.pop();
            for (int v = start; v != -1; v = next[v]) {
                type[v] = curtype;
                vid[v] = pos++;
                inv[vid[v]] = v;
                head[v] = start;
                for (auto e : G[v]) if (e != par[v] && e != next[v]) st.push(e);
            }
        }
    }
    
    // node query [u, v], f([left, right])
    void foreach_nodes(int u, int v, const function<void(int,int)> &f) {
        while (true) {
            if (vid[u] > vid[v]) swap(u, v);
            f(max(vid[head[v]], vid[u]), vid[v]);
            if (head[u] != head[v]) v = par[head[v]];
            else break;
        }
    }
    
    // edge query [u, v], f([left, right])
    void foreach_edges(int u, int v, const function<void(int,int)> &f) {
        while (true) {
            if (vid[u] > vid[v]) swap(u, v);
            if (head[u] != head[v]) {
                f(vid[head[v]], vid[v]);
                v = par[head[v]];
            }
            else {
                if (u != v) {
                    f(vid[u]+1, vid[v]);
                }
                break;
            }
        }
    }

    // LCA
    int get_lca(int u, int v) {
        while (true) {
            if (vid[u] > vid[v]) swap(u, v);
            if (head[u] == head[v]) return u;
            v = par[head[v]];
        }
    }
};



//------------------------------//
// Examples
//------------------------------//

int main () {
    // グラフの入力
    int N;
    cin >> N;
    Graph G(N);
    for (int i = 0; i < N; ++i) {
        int num;
        cin >> num;
        for (int j = 0; j < num; ++j) {
            int c;
            cin >> c;
            G[i].push_back(c);
            G[c].push_back(i);
        }
    }
    
    // HL 分解
    HLDecomposition hld(G);
    hld.build();
    int Q;
    cin >> Q;
    for (int q = 0; q < Q; ++q) {
        int u, v;
        cin >> u >> v;
        cout << hld.get_lca(u, v) << endl;
    }
}
