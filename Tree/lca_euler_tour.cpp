//
// LCA (by Euler Tree)
//
// verified:
//   AOJ Course GRL_5_C Tree - Lowest Common Ancestor
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_5_C&lang=jp
//


#include <bits/stdc++.h>
using namespace std;


// Sparse Table
template<class MeetSemiLattice> struct SparseTable {
    vector<vector<pair<MeetSemiLattice,int> > > dat;
    vector<int> height;
    
    SparseTable() { }
    SparseTable(const vector<MeetSemiLattice> &vec) { init(vec); }
    void init(const vector<MeetSemiLattice> &vec) {
        int n = (int)vec.size(), h = 1;
        while ((1<<h) < n) ++h;
        dat.assign(h, vector<pair<MeetSemiLattice,int> >(1<<h));
        height.assign(n+1, 0);
        for (int i = 2; i <= n; i++) height[i] = height[i>>1]+1;
        for (int i = 0; i < n; ++i) dat[0][i] = {vec[i], i};
        for (int i = 1; i < h; ++i)
            for (int j = 0; j < n; ++j)
                dat[i][j] = min(dat[i-1][j], dat[i-1][min(j+(1<<(i-1)),n-1)]);
    }
    
    pair<MeetSemiLattice,int> get(int a, int b) {
        return min(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};

// Euler Tour
using Graph = vector<vector<int> >;
struct EulerTour {
    // main results
    Graph tree;
    vector<int> depth;
    vector<int> node; // the node-number of i-th element of Euler-tour
    vector<int> vf, ve; // the index of Euler-tour of node v
    vector<int> eid; // the index of edge e (i*2 + (0: dir to leaf, 1: dir to root))
    
    // sub results
    SparseTable<int> st; // depth (to find LCA)
    
    // initialization
    EulerTour(const Graph &tree_) { init(tree_); }
    void init(const Graph &tree_) {
        tree = tree_;
        int V = (int)tree.size();
        depth.resize(V*2-1);
        node.resize(V*2-1);
        vf.resize(V);
        ve.resize(V);
        eid.resize((V-1)*2);
        int k = 0;
        dfs(0, -1, 0, k);
        st.init(depth);
    }
    
    void dfs(int v, int par, int dep, int &ord) {
        node[ord] = v;
        depth[ord] = dep;
        vf[v] = ve[v] = ord;
        ++ord;
        for (auto e : tree[v]) {
            if (e == par) continue;
            dfs(e, v, dep+1, ord);
            node[ord] = v;
            depth[ord] = dep;
            ve[v] = ord;
            ++ord;
        }
    }
    
    inline int get_lca(int u, int v) {
        int a = vf[u], b = vf[v];
        if (a > b) swap(a, b);
        return node[st.get(a, b+1).second];
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
    
    // Euler Tree
    EulerTour et(G);
    int Q;
    cin >> Q;
    for (int q = 0; q < Q; ++q) {
        int u, v;
        cin >> u >> v;
        cout << et.get_lca(u, v) << endl;
    }
}

