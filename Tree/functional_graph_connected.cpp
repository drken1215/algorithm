//
// Functional Graph をサイクルと森に分解する
//
// verified:
//   第二回日本最強プログラマー学生選手権 H - Shipping
//     https://atcoder.jp/contests/jsc2021/tasks/jsc2021_h
//


#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Utility
//------------------------------//

using ll = long long;
using u32 = unsigned int;
using u64 = unsigned long long;
using i128 = __int128_t;
using u128 = __uint128_t;
using pint = pair<int, int>;
using pll = pair<long long, long long>;
using tint = array<int, 3>;
using tll = array<long long, 3>;
using fint = array<int, 4>;
using fll = array<long long, 4>;
using qint = array<int, 5>;
using qll = array<long long, 5>;
using sint = array<int, 6>;
using sll = array<long long, 6>;
using vint = vector<int>;
using vll = vector<long long>;
using dint = deque<int>;
using dll = deque<long long>;
using vvint = vector<vector<int>>;
using vvll = vector<vector<long long>>;
using vpint = vector<pair<int, int>>;
using vpll = vector<pair<long long, long long>>;
template<class T> using min_priority_queue = priority_queue<T, vector<T>, greater<T>>;

template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }
template<class S, class T> inline auto maxll(S a, T b) { return max(ll(a), ll(b)); }
template<class S, class T> inline auto minll(S a, T b) { return min(ll(a), ll(b)); }
template<class T> auto max(const T &a) { return *max_element(a.begin(), a.end()); }
template<class T> auto min(const T &a) { return *min_element(a.begin(), a.end()); }
template<class T> auto argmax(const T &a) { return max_element(a.begin(), a.end()) - a.begin(); }
template<class T> auto argmin(const T &a) { return *min_element(a.begin(), a.end()) - a.begin(); }
template<class T> auto accum(const vector<T> &a) { return accumulate(a.begin(), a.end(), T()); }
template<class T> auto accum(const deque<T> &a) { return accumulate(a.begin(), a.end(), T()); }

#define REP(i, a) for (long long i = 0; i < (long long)(a); i++)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); i++)
#define RREP(i, a) for (long long i = (a)-1; i >= (long long)(0); --i)
#define RREP2(i, a, b) for (long long i = (b)-1; i >= (long long)(a); --i)
#define EB emplace_back
#define PB push_back
#define MP make_pair
#define MT make_tuple
#define FI first
#define SE second
#define ALL(x) x.begin(), x.end()
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl

// input
template<class T> istream& operator >> (istream &is, vector<T> &P)
{ for (int i = 0; i < P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, deque<T> &P)
{ for (int i = 0; i < P.size(); ++i) cin >> P[i]; return is; }

// output
template<class S, class T> ostream& operator << (ostream &s, const pair<S, T> &P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 2> &P)
{ return s << '<' << P[0] << "," << P[1] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 3> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 4> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << "," << P[3] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 5> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << "," << P[3] << "," << P[4] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 6> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << "," << P[3] << "," << P[4] << "," << P[5] << '>'; }
template<class T> ostream& operator << (ostream &s, const vector<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const deque<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const vector<vector<T>> &P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, const set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const multiset<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const unordered_set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const unordered_map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }

// min non-negative i such that n <= 2^i
template<class T> T ceil_pow2(T n) {
    T i = 0;
    while ((T(1) << i) < T(n)) i++;
    return i;
}


//------------------------------//
// Functional Graph
//------------------------------//

// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = -1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
template<class T = long long> struct Graph {
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    
    Graph(int n = 0) : list(n), reversed_list(n) { }
    void init(int n = 0) {
        list.assign(n, vector<Edge<T>>());
        reversed_list.assign(n, vector<Edge<T>>());
    }
    Graph &operator = (const Graph &g) {
        list = g.list, reversed_list = g.reversed_list;
        return *this;
    }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    const vector<Edge<T>> &get_rev_edges(int i) const { return reversed_list[i]; }
    const size_t size() const { return list.size(); }
    const void clear() { list.clear(); }
    const void resize(int n) { list.resize(n); }
        
    void add_edge(int from, int to, T val = 1) {
        list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }
    
    void add_bidirected_edge(int from, int to, T val = 1) {
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
        int D = ceil_pow2(N);
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
    ll N, M; cin >> N >> M;
    Graph<ll> G(N);
    REP(i, N) {
        ll A, C; cin >> A >> C, A--;
        G.add_edge(i, A, C);
    }
    RunConnectedFunctionalGraph<ll> fg(G);
    int C = fg.cycle.size();

    // step 2: ひげ部分の imos 法処理と、サイクル上のパスの取得
    vll imos(N, 0);
    vll X(M), Y(M);
    vector<pll> paths;
    REP(i, M) {
        cin >> X[i] >> Y[i], X[i]--, Y[i]--;
        int rx = fg.root[X[i]], ry = fg.root[Y[i]];
        if (rx == ry) {
            auto lca = fg.get_lca(X[i], Y[i]);
            imos[X[i]]++, imos[lca]--, imos[Y[i]]++, imos[lca]--;
        } else {
            imos[X[i]]++, imos[rx]--, imos[Y[i]]++, imos[ry]--;
            int idx = fg.cmp[rx], idy = fg.cmp[ry];
            if (idx > idy) swap(idx, idy);
            paths.EB(idx, idy);
        }
    }
    auto rec = [&](auto &&rec, int v) -> ll {
        ll res = 0;
        for (auto e : fg.childs[v]) {
            res += rec(rec, e.to);
            imos[v] += imos[e.to];
            if (imos[e.to]) res += e.val;
        }
        return res;
    };
    ll mori = 0;
    for (auto r : fg.roots) mori += rec(rec, r);

    // step 3: サイクル部分の初期コストの計算
    vector<pll> ini(C);
    ll all_sum = 0;
    REP(i, C) ini[i] = pll(0, fg.cycle[i].val), all_sum += fg.cycle[i].val;
    const ll INF = 1LL << 60;
    auto op = [&](pll a, pll b) -> pll {
        if (a.first == b.first) return pll(a.first, a.second + b.second);
        else return min(a, b);
    };
    auto mapping = [&](ll f, pll a) -> pll { return pll(a.first + f, a.second); };
    auto composition = [&](ll g, ll f) -> ll { return g + f; };
    LazySegmentTree<pll, ll> seg(ini, op, mapping, composition, pll(INF, 0), 0);
    auto calc_score = [&]() -> ll {
        auto [mi, sum] = seg.all_prod();
        return (mi == 0 ? all_sum - sum : all_sum);
    };
    vvll left(N), right(N);
    REP(i, paths.size()) {
        auto [x, y] = paths[i];
        left[x].EB(i), right[y].EB(i), seg.apply(x, y, 1);
    }
    ll cur = mori + calc_score();

    // step 4: サイクル部分の差分更新の計算
    ll res = cur;
    REP(i, N-1) {
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
        chmin(res, cur);
    }
    cout << res << endl;
}


int main() {
    The2ndSaikyoKonH();
}