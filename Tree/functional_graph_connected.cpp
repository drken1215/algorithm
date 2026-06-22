//
// Functional Graph をサイクルと森に分解する
//
// verified:
//   第二回日本最強プログラマー学生選手権 H - Shipping
//     https://atcoder.jp/contests/jsc2021/tasks/jsc2021_h
//


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
    int V;
    bool record_reversed_edges = false, record_edge_index = false;
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    vector<unordered_map<int, int>> id;  // id[v][w] := the index of node w in G[v]

    // constructors
    Graph(int n = 0, bool rre = false, bool rei = false) {
        init(n, rre, rei);
    }
    void init(int n = 0, bool rre = false, bool rei = false) {
        V = n, record_reversed_edges = rre, record_edge_index = rei;
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

    // input (only tree-case)
    friend istream& operator >> (istream &is, Graph &G) {
        for (int i = 0; i < G.V - 1; i++) {
            int u, v;
            is >> u >> v, u--, v--;
            G.add_bidirected_edge(u, v);
        }
        return is;
    }

    // output
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

// 連結な Functional Graph を、サイクルと森に分解する
// G[v] の出次数が 1 でなければならない
template<class T = long long> struct RunConnectedFunctionalGraph {
    // cycle
    const int NOT_IN_CYCLE = -1;
    vector<int> roots;  // nodes in the cycle
    vector<Edge<T>> cycle;  // the cycle
    vector<int> cmp;  // order in tye cycle

    // trees
    vector<vector<Edge<T>>> childs;
    vector<unordered_map<int,int>> id;  // id[v][w] := the index of node w in G[v]
    vector<long long> siz;  // the size of v-subtree
    
    // for finding lca
    vector<vector<int>> parent;
    vector<int> root, depth;

    // Euler tour
    vector<int> tour; // the node-number of i-th element of Euler-tour
    vector<int> v_s_id, v_t_id; // the index of Euler-tour of node v
    vector<int> e_id; // the index of edge e (v*2 + (0: root to leaf, 1: leaf to root)

    // constructor
    RunConnectedFunctionalGraph() {}
    RunConnectedFunctionalGraph(const Graph<T> &G, int s = 0) {
        init(G, s);
    }

    // get first / last id of node v in Euler tour
    int vs(int v) { return v_s_id[v]; }
    int vt(int v) { return v_t_id[v]; }
    int get_v(int id) { return tour[id]; }

    // get edge-id of (pv, v) in Euler tour
    int e(int v, bool leaf_to_root = false) {
        assert(cmp[v] == NOT_IN_CYCLE);
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

    // get_parent(v, p) := the parent of v directed for p
    int get_parent(int v) { return parent[0][v];  }
    int get_root(int v) { return root[v]; }
    int kth_ancestor(int v, int k) {
        if (k > depth[v]) return root[v];
        int goal_depth = depth[v] - k;
        for (int i = (int)parent.size()-1; i >= 0; i--)
            if (parent[i][v] != -1 && depth[parent[i][v]] >= goal_depth) 
                v = parent[i][v];
        return v;
    }
    int get_parent(int v, int p) {
        assert(v != p && root[v] == root[p]);
        int lca = get_lca(v, p);
        if (lca != v) return parent[0][v];
        else return kth_ancestor(p, depth[p] - depth[v] - 1);
    }

    // lca(u, v)
    int get_lca(int u, int v) {
        assert(root[u] == root[v]);
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
        if (root[u] == root[v]) {
            int lca = get_lca(u, v);
            return depth[u] + depth[v] - depth[lca] * 2;
        } else {
            int res = depth[u] + depth[v];
            u = root[u], v = root[v];
            int cycledis = max(cmp[u], cmp[v]) - min(cmp[u], cmp[v]);
            cycledis = min(cycledis, (int)cycle.size() - cycledis);
            res += cycledis;
            return res;
        }
    }

    // is node v in s-t path?
    bool is_on_path(int s, int t, int v) {
        return get_dist(s, v) + get_dist(v, t) == get_dist(s, t);
    };
    
    // init
    void init(const Graph<T> &G, int s = 0) {
        int N = (int)G.size();

        // step 0: assertion
        for (int v = 0; v < N; v++) assert(G[v].size() == 1);
        
        // step 1: detect a node in the cycle
        roots.clear(), cycle.clear();
        vector<bool> seen(N, false), finished(N, false);
        int r = s;
        do {
            assert(r != -1);
            seen[r] = true;
            r = G[r][0].to; 
        } while (!seen[r]);

        // step 2: construct cycle
        int v = r, iter = 0;
        cmp.assign(N, NOT_IN_CYCLE);
        do {
            roots.emplace_back(v);
            cycle.emplace_back(G[v][0]);
            cmp[v] = iter++;
            v = G[v][0].to;
        } while (v != r);

        // step 3: construct trees
        int D = 0;
        while ((1LL << D) < N) D++;
        parent.assign(D + 1, vector<int>(N, -1)), childs.resize(N);
        for (int v = 0; v < N; v++) {
            if (cmp[v] != NOT_IN_CYCLE) {
                parent[0][v] = v;
            } else {
                childs[G[v][0].to].emplace_back(Edge<T>(G[v][0].to, v, G[v][0].val));
                parent[0][v] = G[v][0].to;
            }
        }
        for (int i = 0; i < D; i++) for (int v = 0; v < N; v++) {
            parent[i + 1][v] = parent[i][parent[i][v]];
        }

        // step 4: run trees
        depth.resize(N), root.resize(N), siz.resize(N), id.resize(N);
        tour.resize(N * 2 - 1), v_s_id.resize(N), v_t_id.resize(N), e_id.resize(N * 2);
        int ord = 0;
        auto rec = [&](auto &&rec, int v, int d, int r) -> int {
            int sum = 1;
            depth[v] = d, root[v] = r, tour[ord] = v, v_s_id[v] = v_t_id[v] = ord;
            ord++;
            for (int i = 0; i < (int)childs[v].size(); i++) {
                int ch = childs[v][i].to;
                id[v][ch] = i;
                e_id[ch * 2] = ord - 1;
                sum += rec(rec, ch, d + 1, r);
                tour[ord] = v, v_t_id[v] = ord, e_id[ch * 2 + 1] = ord - 1;
                ord++;
            }
            siz[v] = sum;
            return sum;
        };
        for (auto r : roots) rec(rec, r, 0, r);
    }
};


//------------------------------//
// Solver
//------------------------------//


// 第二回日本最強プログラマー学生選手権 H - Shipping
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
    LazySegmentTree(const FuncMonoid op, const FuncAction act, const FuncComposition comp,
                    const Monoid &identity_monoid, const Action &identity_action) 
                    : OP(op), ACT(act), COMP(comp), 
                    IDENTITY_MONOID(identity_monoid), IDENTITY_ACTION(identity_action) {}
    LazySegmentTree(int n, const FuncMonoid op, const FuncAction act, const FuncComposition comp,
                    const Monoid &identity_monoid, const Action &identity_action) {
        init(n, op, act, comp, identity_monoid, identity_action);
    }
    LazySegmentTree(const vector<Monoid> &v,
                    const FuncMonoid op, const FuncAction act, const FuncComposition comp,
                    const Monoid &identity_monoid, const Action &identity_action) {
        init(v, op, act, comp, identity_monoid, identity_action);
    }
    void init(const FuncMonoid op, const FuncAction act, const FuncComposition comp,
              const Monoid &identity_monoid, const Action &identity_action) {
        OP = op, ACT = act, COMP = comp;
        IDENTITY_MONOID = identity_monoid, IDENTITY_ACTION = identity_action;      
    }
    void init(int n) {
        N = n, 
        log = 0, offset = 1;
        while (offset < N) ++log, offset <<= 1;
        dat.assign(offset * 2, IDENTITY_MONOID);
        lazy.assign(offset * 2, IDENTITY_ACTION);
    }
    void init(const vector<Monoid> &v) {
        init((int)v.size());
        build(v);
    }
    void init(int n, const FuncMonoid op, const FuncAction act, const FuncComposition comp,
              const Monoid &identity_monoid, const Action &identity_action) {
        init(op, act, comp, identity_monoid, identity_action);
        init(n);
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

void The2ndSaikyoKonH() {
    // step 1: 入力と、Functional Graph グラフの分解（無向グラフだが、Functional Graph として見た方がやりやすい）
    long long N, M; 
    cin >> N >> M;
    Graph<long long> G(N);
    for (int i = 0; i < N; i++) {
        long long A, C; 
        cin >> A >> C, A--;
        G.add_edge(i, A, C);
    }
    RunConnectedFunctionalGraph<long long> fg(G);
    int C = fg.cycle.size();

    // step 2: ひげ部分の imos 法処理と、サイクル上のパスの取得
    vector<long long> imos(N, 0), X(M), Y(M);
    vector<pair<long long, long long>> paths;
    for (int i = 0; i < M; i++) {
        cin >> X[i] >> Y[i], X[i]--, Y[i]--;
        int rx = fg.root[X[i]], ry = fg.root[Y[i]];
        if (rx == ry) {
            auto lca = fg.get_lca(X[i], Y[i]);
            imos[X[i]]++, imos[lca]--, imos[Y[i]]++, imos[lca]--;
        } else {
            imos[X[i]]++, imos[rx]--, imos[Y[i]]++, imos[ry]--;
            int idx = fg.cmp[rx], idy = fg.cmp[ry];
            if (idx > idy) swap(idx, idy);
            paths.emplace_back(idx, idy);
        }
    }
    auto rec = [&](auto &&rec, int v) -> long long {
        long long res = 0;
        for (auto e : fg.childs[v]) {
            res += rec(rec, e.to);
            imos[v] += imos[e.to];
            if (imos[e.to]) res += e.val;
        }
        return res;
    };
    long long mori = 0;
    for (auto r : fg.roots) mori += rec(rec, r);

    // step 3: サイクル部分の初期コストの計算
    using pll = pair<long long, long long>;
    vector<pair<long long, long long>> ini(C);
    long long all_sum = 0;
    for (int i = 0; i < C; i++) {
        ini[i] = make_pair(0, fg.cycle[i].val); 
        all_sum += fg.cycle[i].val;
    }
    const long long INF = 1LL << 60;
    auto op = [&](pll a, pll b) -> pll {
        if (a.first == b.first) return pll(a.first, a.second + b.second);
        else return min(a, b);
    };
    auto mapping = [&](long long f, pll a) -> pll { return pll(a.first + f, a.second); };
    auto composition = [&](long long g, long long f) -> long long { return g + f; };
    LazySegmentTree<pll, long long> seg(ini, op, mapping, composition, pll(INF, 0), 0);
    auto calc_score = [&]() -> long long {
        auto [mi, sum] = seg.all_prod();
        return (mi == 0 ? all_sum - sum : all_sum);
    };
    vector<vector<long long>> left(N), right(N);
    for (int i = 0; i < paths.size(); i++) {
        auto [x, y] = paths[i];
        left[x].emplace_back(i), right[y].emplace_back(i), seg.apply(x, y, 1);
    }
    long long cur = mori + calc_score();

    // step 4: サイクル部分の差分更新の計算
    long long res = cur;
    for (int i = 0; i < N - 1; i++) {
        // erase を ((i-1+N)%N, i) から (i, i+1) に移す
        // 始点が i になっているやつは、反転する
        for (auto id : left[i]) {
            auto [x, y] = paths[id];  // (x < y)
            seg.apply(x, y, -1);
            seg.apply(y, C, 1), seg.apply(0, x, 1);
        }
        // 終点が i になっているやつは、反転する
        for (auto id : right[i]) {
            auto [x, y] = paths[id];  // (x < y)
            seg.apply(x, y, 1);
            seg.apply(y, C, -1), seg.apply(0, x, -1);
        }
        cur = mori + calc_score();
        res = min(res, cur);
    }
    cout << res << endl;
}


int main() {
    The2ndSaikyoKonH();
}