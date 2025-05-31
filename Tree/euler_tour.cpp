//
// Euler Tour (update for nodes)
//
// verified:
//   ABC 406 F - Compare Tree Weights
//     https://atcoder.jp/contests/abc406/tasks/abc406_f
//
//   ABC 294 G - Distance Queries on a Tree
//     https://atcoder.jp/contests/abc294/tasks/abc294_g
//
//   AOJ 2667 Tree
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2667
//
//   ABC 133 F - Colorful Tree
//     https://atcoder.jp/contests/abc133/tasks/abc133_f
//

/*
    ・頂点の行きがけ順の取得：vs(v)
    ・頂点の帰りがけ順の取得：vt(v)
    ・辺 (p(v), v) の取得：e(v, false)（p(v) は v の親)
    ・辺 (v, p(v)) の取得：e(v, true)（p(v) は v の親)
    ・パス 0-v クエリ：区間 [0, vs(v)) への処理
    ・v-部分木 クエリ：区間 [vs(v), vt(v) + 1) への処理
*/


#include <bits/stdc++.h>
using namespace std;


// Run Tree (including Euler Tour)
template<class Graph = vector<vector<int>>> struct RunTree {
    // id[v][w] := the index of node w in G[v]
    vector<unordered_map<int, int>> id;

    // num[v][i] := the size of subtree of G[v][i] with parent v
    vector<vector<long long>> num;
    
    // for finding lca
    int root;
    vector<vector<int>> parent;
    vector<int> depth;

    // Euler tour
    vector<int> tour; // the node-number of i-th element of Euler-tour
    vector<int> v_s_id, v_t_id; // the index of Euler-tour of node v
    vector<int> e_id; // the index of edge e (v*2 + (0: root to leaf, 1: leaf to root))

    // constructor
    RunTree() {}
    RunTree(const Graph &G, int root = 0) : root(root) {
        init(G, root);
    }
    
    // init
    void init(const Graph &G, int root = 0) {
        int N = (int)G.size();
        id.assign(N, unordered_map<int,int>()), num.assign(N, vector<long long>());
        for (int v = 0; v < N; v++) num[v].assign((int)G[v].size(), 0);
        int h = 1, ord = 0;
        while ((1<<h) < N) h++;
        parent.assign(h, vector<int>(N, -1)), depth.resize(N);
        tour.resize(N*2-1), v_s_id.resize(N), v_t_id.resize(N), e_id.resize(N*2);
        rec(G, root, -1, 0, ord);
        for (int i = 0; i+1 < (int)parent.size(); ++i) {
            for (int v = 0; v < N; v++)
                if (parent[i][v] != -1)
                    parent[i+1][v] = parent[i][parent[i][v]];
        }
    }

    // get_size(u, v) := the size of subtree v with parent u
    long long get_size(int u, int v) {
        return num[u][id[u][v]];
    }

    // get first / last id of node v in Euler tour
    int vs(int v) { return v_s_id[v]; }
    int vt(int v) { return v_t_id[v]; }
    int get_v(int id) { return tour[id]; }

    // get edge-id of (pv, v) in Euler tour
    int e(int v, bool leaf_to_root = false) {
        assert(v != root);
        if (!leaf_to_root) return e_id[v * 2];
        else return e_id[v * 2 + 1];
    }
    int e(int u, int v) {
        if (depth[u] < depth[v]) return e(v);
        else return e(u, false);
    }
    pair<int, int> get_e(int id) { 
        return make_pair(tour[id], tour[id + 1]);
    }

    // lca(u, v)
    int get_lca(int u, int v) {
        if (depth[u] > depth[v]) swap(u, v);
        for (int i = 0; i < (int)parent.size(); i++) {
            if ((depth[v] - depth[u]) & (1<<i))
                v = parent[i][v];
        }
        if (u == v) return u;
        for (int i = (int)parent.size()-1; i >= 0; i--) {
            if (parent[i][u] != parent[i][v]) {
                u = parent[i][u];
                v = parent[i][v];
            }
        }
        return parent[0][u];
    }

    // dist(u, v)
    long long get_dist(int u, int v) {
        int lca = get_lca(u, v);
        return depth[u] + depth[v] - depth[lca]*2;
    }

    // get_parent(v, p) := the parent of v directed for p
    int get_parent(int v, int p) {
        if (v == p) return -1;
        int lca = get_lca(v, p);
        if (lca != v) return parent[0][v];
        for (int i = (int)parent.size()-1; i >= 0; i--) {
            if (parent[i][p] != -1 && depth[parent[i][p]] > depth[v]) {
                p = parent[i][p];
            }
        }
        return p;
    }
    
    // rec
    int rec(const Graph &G, int v, int p, int d, int &ord) {
        int p_index = -1;
        int sum = 1;
        parent[0][v] = p, depth[v] = d;
        tour[ord] = v, v_s_id[v] = v_t_id[v] = ord;
        ord++;
        for (int i = 0; i < (int)G[v].size(); i++) {
            int ch = G[v][i];
            id[v][ch] = i;
            if (ch == p) {
                p_index = i;
                continue;
            }
            e_id[ch * 2] = ord - 1;
            int s = rec(G, ch, v, d+1, ord);
            num[v][i] = s;
            sum += s;
            tour[ord] = v;
            v_t_id[v] = ord;
            e_id[ch * 2 + 1] = ord - 1;
            ord++;
        }
        if (p_index != -1) num[v][p_index] = (int)G.size() - sum;
        return sum;
    }
};


//------------------------------//
// Examples
//------------------------------//

// ABC 406 F - Compare Tree Weights
template <class Abel> struct BIT {
    Abel UNITY_SUM = 0;
    vector<Abel> dat;
    
    // [0, n)
    BIT(int n, Abel unity = 0) : UNITY_SUM(unity), dat(n, unity) { }
    void init(int n) {
        dat.assign(n, UNITY_SUM);
    }
    int size() const {
        return (int)dat.size();
    }
    
    // a is 0-indexed
    inline void add(int a, Abel x) {
        for (int i = a; i < (int)dat.size(); i |= i + 1)
            dat[i] = dat[i] + x;
    }
    
    // [0, a), a is 0-indexed, [a, b), a and b are 0-indexed
    inline Abel sum(int a) const {
        Abel res = UNITY_SUM;
        for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[i];
        return res;
    }
    inline Abel sum(int a, int b) const {
        return sum(b) - sum(a);
    }
    inline Abel operator [] (int i) const {
        return sum(i, i + 1);
    }
    
    // debug
    friend ostream& operator << (ostream &s, const BIT &bit) {
        for (int i = 0; i < (int)bit.size(); ++i) s << bit[i] << " ";
        return s;
    }
};

void ABC_406_F() {
    using Graph = vector<vector<int>>;
    int N, Q, u, v;
    cin >> N;
    Graph G(N);
    vector<pair<int,int>> edges(N-1);
    for (int i = 0; i < N-1; i++) {
        cin >> u >> v, u--, v--;
        edges[i] = {u, v};
        G[u].push_back(v), G[v].push_back(u);
    }
    RunTree rt(G);
    cin >> Q;
    long long all = 0;
    BIT<long long> bit(rt.tour.size() + 1);
    for (int qid = 0; qid < Q; qid++) {
        long long type, v, x, y;
        cin >> type;
        if (type == 1) {
            cin >> v >> x, v--;
            bit.add(rt.vs(v), x);
            all += x;
        } else if (type == 2) {
            cin >> y, y--;
            int u = edges[y].first, v = edges[y].second;
            if (rt.depth[u] > rt.depth[v]) swap(u, v);

            long long uv_size = rt.get_size(u, v), vu_size = rt.get_size(v, u);
            long long uv_sum = bit.sum(rt.vs(v), rt.vt(v)+1);  // +1 is necessary!
            long long vu_sum = all - uv_sum;
            long long uv = uv_size + uv_sum, vu = vu_size + vu_sum;
            long long res = abs(uv - vu);
            cout << res << '\n';
        } 
    }
}


// ABC 294 G - Distance Queries on a Tree
void ABC_294_G() {
    long long N, Q, typ;
    cin >> N;
    vector<vector<int>> G(N);
    vector<array<long long, 3>> edges(N-1);
    for (int i = 0; i < N-1; i++) {
        long long u, v, w;
        cin >> u >> v >> w, u--, v--;
        G[u].emplace_back(v), G[v].emplace_back(u);
        edges[i] = array<long long, 3>({u, v, w});
    }
    RunTree rt(G);
    BIT<long long> bit(N * 2);
    for (auto [u, v, w] : edges) {
        if (rt.depth[u] > rt.depth[v]) swap(u, v);
        bit.add(rt.e(v, false), w);
        bit.add(rt.e(v, true), -w);
    }

    cin >> Q;
    while (Q--) {
        cin >> typ;
        if (typ == 1) {
            long long i, w;
            cin >> i >> w, i--;
            auto [u, v, pw] = edges[i];
            if (rt.depth[u] > rt.depth[v]) swap(u, v);
            int e1 = rt.e(v, false), e2 = rt.e(v, true);
            bit.add(e1, w - bit[e1]);
            bit.add(e2, -w - bit[e2]);
        } else {
            int u, v;
            cin >> u >> v, u--, v--;
            int l = rt.get_lca(u, v);
            long long res = bit.sum(0, rt.vs(u)) + bit.sum(0, rt.vs(v)) 
                - bit.sum(0, rt.vs(l)) * 2;
            cout << res << '\n';
        }
    }
}


// AOJ 2667 Tree
// Lazy Segment Tree
template<class Monoid, class Action> struct LazySegmentTree {
    // various function types
    using FuncMonoid = function<Monoid(Monoid, Monoid)>;
    using FuncAction = function<Monoid(Action, Monoid)>;
    using FuncComposition = function<Action(Action, Action)>;

    // core member
    int N;
    FuncMonoid OP;
    FuncAction ACT;
    FuncComposition COMP;
    Monoid IDENTITY_MONOID;
    Action IDENTITY_ACTION;
    
    // inner data
    int log, offset;
    vector<Monoid> dat;
    vector<Action> lazy;
    
    // constructor
    LazySegmentTree() {}
    LazySegmentTree(int n, const FuncMonoid op, const FuncAction act, const FuncComposition comp,
                    const Monoid &identity_monoid, const Action &identity_action) {
        init(n, op, act, comp, identity_monoid, identity_action);
    }
    LazySegmentTree(const vector<Monoid> &v,
                    const FuncMonoid op, const FuncAction act, const FuncComposition comp,
                    const Monoid &identity_monoid, const Action &identity_action) {
        init(v, op, act, comp, identity_monoid, identity_action);
    }
    void init(int n, const FuncMonoid op, const FuncAction act, const FuncComposition comp,
              const Monoid &identity_monoid, const Action &identity_action) {
        N = n, OP = op, ACT = act, COMP = comp;
        IDENTITY_MONOID = identity_monoid, IDENTITY_ACTION = identity_action;
        log = 0, offset = 1;
        while (offset < N) ++log, offset <<= 1;
        dat.assign(offset * 2, IDENTITY_MONOID);
        lazy.assign(offset * 2, IDENTITY_ACTION);
    }
    void init(const vector<Monoid> &v,
              const FuncMonoid op, const FuncAction act, const FuncComposition comp,
              const Monoid &identity_monoid, const Action &identity_action) {
        init((int)v.size(), op, act, comp, identity_monoid, identity_action);
        build(v);
    }
    void build(const vector<Monoid> &v) {
        assert(N == (int)v.size());
        for (int i = 0; i < N; ++i) dat[i + offset] = v[i];
        for (int k = offset - 1; k > 0; --k) pull_dat(k);
    }
    int size() const {
        return N;
    }
    
    // basic functions for lazy segment tree
    void pull_dat(int k) {
        dat[k] = OP(dat[k * 2], dat[k * 2 + 1]);
    }
    void apply_lazy(int k, const Action &f) {
        dat[k] = ACT(f, dat[k]);
        if (k < offset) lazy[k] = COMP(f, lazy[k]);
    }
    void push_lazy(int k) {
        apply_lazy(k * 2, lazy[k]);
        apply_lazy(k * 2 + 1, lazy[k]);
        lazy[k] = IDENTITY_ACTION;
    }
    void pull_dat_deep(int k) {
        for (int h = 1; h <= log; ++h) pull_dat(k >> h);
    }
    void push_lazy_deep(int k) {
        for (int h = log; h >= 1; --h) push_lazy(k >> h);
    }
    
    // setter and getter, update A[i], i is 0-indexed, O(log N)
    void set(int i, const Monoid &v) {
        assert(0 <= i && i < N);
        int k = i + offset;
        push_lazy_deep(k);
        dat[k] = v;
        pull_dat_deep(k);
    }
    Monoid get(int i) {
        assert(0 <= i && i < N);
        int k = i + offset;
        push_lazy_deep(k);
        return dat[k];
    }
    Monoid operator [] (int i) {
        return get(i);
    }
    
    // apply f for index i
    void apply(int i, const Action &f) {
        assert(0 <= i && i < N);
        int k = i + offset;
        push_lazy_deep(k);
        dat[k] = ACT(f, dat[k]);
        pull_dat_deep(k);
    }
    // apply f for interval [l, r)
    void apply(int l, int r, const Action &f) {
        assert(0 <= l && l <= r && r <= N);
        if (l == r) return;
        l += offset, r += offset;
        for (int h = log; h >= 1; --h) {
            if (((l >> h) << h) != l) push_lazy(l >> h);
            if (((r >> h) << h) != r) push_lazy((r - 1) >> h);
        }
        int original_l = l, original_r = r;
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) apply_lazy(l++, f);
            if (r & 1) apply_lazy(--r, f);
        }
        l = original_l, r = original_r;
        for (int h = 1; h <= log; ++h) {
            if (((l >> h) << h) != l) pull_dat(l >> h);
            if (((r >> h) << h) != r) pull_dat((r - 1) >> h);
        }
    }
    
    // get prod of interval [l, r)
    Monoid prod(int l, int r) {
        assert(0 <= l && l <= r && r <= N);
        if (l == r) return IDENTITY_MONOID;
        l += offset, r += offset;
        for (int h = log; h >= 1; --h) {
            if (((l >> h) << h) != l) push_lazy(l >> h);
            if (((r >> h) << h) != r) push_lazy(r >> h);
        }
        Monoid val_left = IDENTITY_MONOID, val_right = IDENTITY_MONOID;
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) val_left = OP(val_left, dat[l++]);
            if (r & 1) val_right = OP(dat[--r], val_right);
        }
        return OP(val_left, val_right);
    }
    Monoid all_prod() {
        return dat[1];
    }
    
    // get max r that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int max_right(const function<bool(Monoid)> f, int l = 0) {
        if (l == N) return N;
        l += offset;
        push_lazy_deep(l);
        Monoid sum = IDENTITY_MONOID;
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(OP(sum, dat[l]))) {
                while (l < offset) {
                    push_lazy(l);
                    l = l * 2;
                    if (f(OP(sum, dat[l]))) {
                        sum = OP(sum, dat[l]);
                        ++l;
                    }
                }
                return l - offset;
            }
            sum = OP(sum, dat[l]);
            ++l;
        } while ((l & -l) != l);  // stop if l = 2^e
        return N;
    }

    // get min l that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int min_left(const function<bool(Monoid)> f, int r = -1) {
        if (r == 0) return 0;
        if (r == -1) r = N;
        r += offset;
        push_lazy_deep(r - 1);
        Monoid sum = IDENTITY_MONOID;
        do {
            --r;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(OP(dat[r], sum))) {
                while (r < offset) {
                    push_lazy(r);
                    r = r * 2 + 1;
                    if (f(OP(dat[r], sum))) {
                        sum = OP(dat[r], sum);
                        --r;
                    }
                }
                return r + 1 - offset;
            }
            sum = OP(dat[r], sum);
        } while ((r & -r) != r);
        return 0;
    }
    
    // debug stream
    friend ostream& operator << (ostream &s, LazySegmentTree seg) {
        for (int i = 0; i < (int)seg.size(); ++i) {
            s << seg[i];
            if (i != (int)seg.size() - 1) s << " ";
        }
        return s;
    }
    
    // dump
    void dump() {
        for (int i = 0; i <= log; ++i) {
            for (int j = (1 << i); j < (1 << (i + 1)); ++j) {
                cout << "{" << dat[j] << "," << lazy[j] << "} ";
            }
            cout << endl;
        }
    }
};

void AOJ_2667() {
    int N, Q;
    cin >> N >> Q;
    vector<vector<int>> G(N);
    for (int i = 0; i < N-1; ++i) {
        int a, b;
        cin >> a >> b;
        G[a].push_back(b);
        G[b].push_back(a);
    }

    RunTree rt(G);
    using Node = pair<long long, int>;
    auto fm = [&](Node a, Node b) { return Node(a.first + b.first, a.second + b.second); };
    auto fa = [&](long long d, Node a) { a.first += d * a.second; return a; };
    auto fl = [&](long long d, long long e) { return d + e; };
    LazySegmentTree<Node, long long> seg(N*2, fm, fa, fl, Node(0, 0), 0);
    for (int v = 1; v < N; v++) {
        seg.set(rt.e(v, false), Node(0, 1));
        seg.set(rt.e(v, true), Node(0, -1));
    }
    
    for (int q = 0; q < Q; ++q) {
        int type;
        cin >> type;
        if (type == 0) {
            int u, v;
            cin >> u >> v;
            int l = rt.get_lca(u, v);
            long long res = seg.prod(0, rt.vs(u)).first
             + seg.prod(0, rt.vs(v)).first
             - seg.prod(0, rt.vs(l)).first * 2;
            cout << res << '\n';
        } else {
            int u, x;
            cin >> u >> x;
            seg.apply(rt.vs(u), rt.vt(u), x);
        }
    }
}


// ABC 133 F - Colorful Tree
void ABC_133_F() {
    using pint = pair<int, int>;
    using fll = array<long long, 4>;
    int N, Q, a, b, c, d, u, v, w;
    cin >> N >> Q;
    vector<vector<int>> G(N);
    vector<fll> edges(N-1);
    for (int i = 0; i < N-1; i++) {
        cin >> a >> b >> c >> d, a--, b--, c--;
        G[a].emplace_back(b), G[b].emplace_back(a);
        edges[i] = fll({a, b, c, d});
    }
    RunTree rt(G);
    vector<long long> col(N * 2 + 1), num(N * 2 + 1), val(N * 2 + 1);
    for (auto [u, v, c, d] : edges) {
        if (rt.depth[u] > rt.depth[v]) swap(u, v);
        int e1 = rt.e(v, false), e2 = rt.e(v, true);
        col[e1] = col[e2] = c, num[e1] = 1, num[e2] = -1, val[e1] = d, val[e2] = -d;
    }
    vector<vector<fll>> qs(N * 2 + 1);
    for (int qid = 0; qid < Q; qid++) {
        cin >> c >> w >> u >> v, c--, u--, v--;
        long long l = rt.get_lca(u, v);
        qs[rt.vs(u)].emplace_back(fll({c, w, 1, qid}));
        qs[rt.vs(v)].emplace_back(fll({c, w, 1, qid}));
        qs[rt.vs(l)].emplace_back(fll({c, w, -2, qid}));
    }
    long long sum = 0;
    vector<long long> res(Q, 0), cnum(N+1, 0), csum(N+1, 0);
    for (int id = 0; id < N * 2; id++) {
        sum += val[id], cnum[col[id]] += num[id], csum[col[id]] += val[id];
        for (auto [c, w, factor, qid] : qs[id+1]) {
            res[qid] += (sum - csum[c] + cnum[c] * w) * factor;
        }
    }
    for (int qid = 0; qid < Q; qid++) cout << res[qid] << '\n';
}


int main () {
    //ABC_406_F();
    //ABC_294_G();
    //AOJ_2667();
    ABC_133_F();
}