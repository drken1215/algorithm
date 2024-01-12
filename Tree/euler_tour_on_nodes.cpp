//
// Euler Tour (update for nodes)
//
// verified:
//   AOJ 2667 Tree
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2667
//

#include <iostream>
#include <functional>
#include <vector>
#include <queue>
#include <stack>
using namespace std;


// Sparse Table
template<class MeetSemiLattice> struct SparseTable {
    vector<vector<pair<MeetSemiLattice,int> > > dat;
    vector<int> height;
    
    SparseTable() { }
    SparseTable(const vector<MeetSemiLattice> &vec) { init(vec); }
    void init(const vector<MeetSemiLattice> &vec) {
        int n = (int)vec.size(), h = 0;
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

// Segment Tree
template<class Monoid, class Action> struct SegTree {
    using FuncMonoid = function< Monoid(Monoid, Monoid) >;
    using FuncAction = function< void(Monoid&, Action) >;
    using FuncLazy = function< void(Action&, Action) >;
    FuncMonoid FM;
    FuncAction FA;
    FuncLazy FL;
    Monoid UNITY_MONOID;
    Action UNITY_LAZY;
    int SIZE, HEIGHT;
    vector<Monoid> dat;
    vector<Action> lazy;
    
    SegTree() { }
    SegTree(int n, const FuncMonoid fm, const FuncAction fa, const FuncLazy fl,
            const Monoid &unity_monoid, const Action &unity_lazy)
    : FM(fm), FA(fa), FL(fl), UNITY_MONOID(unity_monoid), UNITY_LAZY(unity_lazy) {
        SIZE = 1; HEIGHT = 0;
        while (SIZE < n) SIZE <<= 1, ++HEIGHT;
        dat.assign(SIZE * 2, UNITY_MONOID);
        lazy.assign(SIZE * 2, UNITY_LAZY);
    }
    void init(int n, const FuncMonoid fm, const FuncAction fa, const FuncLazy fl,
              const Monoid &unity_monoid, const Action &unity_lazy) {
        FM = fm; FA = fa; FL = fl;
        UNITY_MONOID = unity_monoid; UNITY_LAZY = unity_lazy;
        SIZE = 1; HEIGHT = 0;
        while (SIZE < n) SIZE <<= 1, ++HEIGHT;
        dat.assign(SIZE * 2, UNITY_MONOID);
        lazy.assign(SIZE * 2, UNITY_LAZY);
    }
    
    /* set, a is 0-indexed */
    void set(int a, const Monoid &v) { dat[a + SIZE] = v; }
    void build() {
        for (int k = SIZE - 1; k > 0; --k)
            dat[k] = FM(dat[k*2], dat[k*2+1]);
    }
    
    /* update [a, b) */
    inline void evaluate(int k) {
        if (lazy[k] == UNITY_LAZY) return;
        if (k < SIZE) FL(lazy[k*2], lazy[k]), FL(lazy[k*2+1], lazy[k]);
        FA(dat[k], lazy[k]);
        lazy[k] = UNITY_LAZY;
    }
    inline void update(int a, int b, const Action &v, int k, int l, int r) {
        evaluate(k);
        if (a <= l && r <= b)  FL(lazy[k], v), evaluate(k);
        else if (a < r && l < b) {
            update(a, b, v, k*2, l, (l+r)>>1), update(a, b, v, k*2+1, (l+r)>>1, r);
            dat[k] = FM(dat[k*2], dat[k*2+1]);
        }
    }
    inline void update(int a, int b, const Action &v) { update(a, b, v, 1, 0, SIZE); }
    
    /* get [a, b) */
    inline Monoid get(int a, int b, int k, int l, int r) {
        evaluate(k);
        if (a <= l && r <= b)
            return dat[k];
        else if (a < r && l < b)
            return FM(get(a, b, k*2, l, (l+r)>>1), get(a, b, k*2+1, (l+r)>>1, r));
        else
            return UNITY_MONOID;
    }
    inline Monoid get(int a, int b) { return get(a, b, 1, 0, SIZE); }
    inline Monoid operator [] (int a) { return get(a, a+1); }
    
    /* debug */
    void print() {
        for (int i = 0; i < SIZE; ++i) { cout << (*this)[i]; if (i != SIZE) cout << ","; }
        cout << endl;
    }
};

// Euler Tour
using Graph = vector<vector<int> >;
using Node = pair<long long, int>;
const auto fm = [](Node a, Node b) { return Node(a.first + b.first, a.second + b.second); };
const auto fa = [](Node &a, long long d) { a.first += d * a.second; };
const auto fl = [](long long &d, long long e) { d += e; };

struct EulerTour {
    // main results
    Graph tree;
    vector<int> depth;
    vector<int> node; // the node-number of i-th element of Euler-tour
    vector<int> vf, ve; // the index of Euler-tour of node v
    vector<int> eid; // the index of edge e (i*2 + (0: dir to leaf, 1: dir to root))
    
    // sub results
    SparseTable<int> st; // depth (to find LCA)
    
    // segtree
    SegTree<Node, long long> seg;
    
    // initialization
    EulerTour(const Graph &tree_) { init(tree_); }
    void init(const Graph &tree_) {
        tree = tree_;
        int V = (int)tree.size();
        depth.resize(V*2-1); node.resize(V*2-1); vf.resize(V); ve.resize(V); eid.resize((V-1)*2);
        seg.init((V-1)*2, fm, fa, fl, Node(0, 0), 0);
        int k = 0;
        dfs(0, -1, 0, k);
        st.init(depth);
        seg.build();
    }
    
    void dfs(int v, int par, int dep, int &ord) {
        node[ord] = v;
        depth[ord] = dep;
        vf[v] = ve[v] = ord;
        ++ord;
        for (auto e : tree[v]) {
            if (e == par) continue;
            seg.set(ord-1, Node(0, 1));
            dfs(e, v, dep+1, ord);
            node[ord] = v;
            depth[ord] = dep;
            ve[v] = ord;
            seg.set(ord-1, Node(0, -1));
            ++ord;
        }
    }
    
    inline int LCA(int u, int v) {
        int a = vf[u], b = vf[v];
        if (a > b) swap(a, b);
        return node[st.get(a, b+1).second];
    }
    
    inline void update(int v, long long x) {
        seg.update(vf[v], ve[v], x);
    }
    
    inline long long get(int v) {
        return seg.get(0, vf[v]).first;
    }
    
    inline long long get(int u, int v) {
        int lca = LCA(u, v);
        return get(u) + get(v) - get(lca)*2;
    }
};



//------------------------------//
// Examples
//------------------------------//

int main () {
    int N, Q;
    cin >> N >> Q;
    Graph G(N);
    for (int i = 0; i < N-1; ++i) {
        int a, b; scanf("%d %d", &a, &b);
        G[a].push_back(b);
        G[b].push_back(a);
    }

    EulerTour et(G);
    
    for (int q = 0; q < Q; ++q) {
        int type; cin >> type;
        if (type == 0) {
            int u, v; scanf("%d %d", &u, &v);
            long long res = et.get(u, v);
            printf("%lld\n", res);
        }
        else {
            int u, x; scanf("%d %d", &u, &x);
            et.update(u, x);
        }
    }
}
