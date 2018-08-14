//
// Heavy-Light Decomposition
//
// verified:
//   AOJ 2667 Tree
//

#include <iostream>
#include <functional>
#include <vector>
#include <queue>
#include <stack>
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

typedef vector<vector<int> > Graph;
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
    int lca(int u, int v) {
        while (true) {
            if (vid[u] > vid[v]) swap(u, v);
            if (head[u] == head[v]) return u;
            v = par[head[v]];
        }
    }
};

// BIT
template <class Abel> struct BIT {
    vector<Abel> dat[2];
    Abel UNITY_SUM = 0;						// to be set
    
    /* [1, n] */
    BIT(int n) { init(n); }
    void init(int n) { for (int iter = 0; iter < 2; ++iter) dat[iter].assign(n + 1, UNITY_SUM); }

    /* a is 1-indexed */
    inline void sub_add(int p, int a, Abel x) {
        for (int i = a; i < (int)dat[p].size(); i += i & -i)
            dat[p][i] = dat[p][i] + x;
    }
    inline void add(int a, int b, Abel x) {
        sub_add(0, a, x * -(a - 1)); sub_add(1, a, x); sub_add(0, b, x * (b - 1)); sub_add(1, b, x * (-1));
    }
    
    /* [1, a], a is 1-indexed */
    inline Abel sub_sum(int p, int a) {
        Abel res = UNITY_SUM;
        for (int i = a; i > 0; i -= i & -i) res = res + dat[p][i];
        return res;
    }
    inline Abel sum(int a, int b) {
        return sub_sum(0, b - 1) + sub_sum(1, b - 1) * (b - 1) - sub_sum(0, a - 1) - sub_sum(1, a - 1) * (a - 1);
    }
    
    /* debug */
    void print() {
        for (int i = 1; i < (int)dat[0].size(); ++i) cout << sum(i, i + 1) << ",";
        cout << endl;
    }
};

int main() {
    int N, Q; cin >> N >> Q;
    Graph G(N);
    for (int i = 0; i < N-1; ++i) {
        int u, v; scanf("%d %d", &u, &v);
        G[u].push_back(v); G[v].push_back(u);
    }
    HLDecomposition hld(G);
    hld.build();
    BIT<long long> bit(N+1);
    for (int q = 0; q < Q; ++q) {
        int type, a, b; scanf("%d %d %d", &type, &a, &b);
        if (type == 0) {
            long long res = 0;

			// 区間和取得
            hld.foreach_edges(a, b, [&](int l, int r){ res += bit.sum(l, r+1); });
            cout << res << endl;
        }
        else {

		    // 部分木のなす区間に b を加算
            bit.add(hld.vid[a]+1, hld.vid[a]+hld.subsize[a], b);
        }
    }
}
