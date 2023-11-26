//
// Lazy Segment Tree
//
// verified:
//   AOJ Course Range Query - RMQ and RUQ (change the val / min)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_F&lang=ja
//
//   AOJ Course Range Query - RSQ and RAQ (add / sum)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_G&lang=ja
//
//   AOJ Course Range Query - RMQ and RAQ (add / min) - Starry Sky Tree
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_H&lang=ja
//
//   AOJ Course Range Query - RSQ and RUQ (change the val / sum)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_I&lang=ja
//
//   AtCoder ACL Beginner Contest E - Replace Digits
//     https://atcoder.jp/contests/abl/tasks/abl_e
//
//   AtCoder ABC 322 F - Vacation Query
//     https://atcoder.jp/contests/abc322/tasks/abc322_f
//
//   AtCoder ABC 330 E - Mex and Update (for max_right)
//     https://atcoder.jp/contests/abc330/tasks/abc330_e
//


/*
    Lazy Segment Tree は「作用つきモノイド」上で定義される
        OP(x, y): 2 つのモノイド間に定義される演算
        MAPPING(f, x): f(x), モノイド元 x への作用素 f による作用
        COMPOSITION(g, f): g o f, 作用素 f への作用素 g の合成
        IDENTITY_MONOID: モノイドの単位元
        IDENTITY_ACTION: 作用素の単位元

    // Construction
    LazySegmentTree(N, op, mapping, composition, identity_monoid, identity_lazy)
    LazySegmentTree(v, op, mapping, composition, identity_monoid, identity_lazy)
      ex: starry sky tree (区間加算、区間min取得)
        auto op = [](long long x, long long y) { return min(x, y); };
        auto mapping = [](int f, long long x) { f + x; };  // replace x with f + x
        auto composition = [](int g, int f) { g + f; };  // replace f with g + f
        LazySegmentTree<long long, long long> seg(N, op, mapping, composition, (1LL<<60), 0);

    // Queries
    set(i, v): 添字 i の箇所を値 v にする, O(log N)
    apply(l, r, a): 区間 [l, r) を作用素 a を用いて更新する, O(log N)
    prod(l, r): 区間 [l, r) についての演算結果を返す, O(log N)
*/


#include <bits/stdc++.h>
using namespace std;


// Lazy Segment Tree
template<class Monoid, class Action> struct LazySegmentTree {
    // various function types
    using FuncOperator = function<Monoid(Monoid, Monoid)>;
    using FuncMapping = function<Monoid(Action, Monoid)>;
    using FuncComposition = function<Action(Action, Action)>;

    // core member
    int N;
    FuncOperator OP;
    FuncMapping MAPPING;
    FuncComposition COMPOSITION;
    Monoid IDENTITY_MONOID;
    Action IDENTITY_ACTION;
    
    // inner data
    int log, offset;
    vector<Monoid> dat;
    vector<Action> lazy;
    
    // constructor
    LazySegmentTree() {}
    LazySegmentTree(int n,
                    const FuncOperator op,
                    const FuncMapping mapping,
                    const FuncComposition composition,
                    const Monoid &identity_monoid,
                    const Action &identity_action) {
        init(n, op, mapping, composition, identity_monoid, identity_action);
    }
    LazySegmentTree(const vector<Monoid> &v,
                    const FuncOperator op,
                    const FuncMapping mapping,
                    const FuncComposition composition,
                    const Monoid &identity_monoid,
                    const Action &identity_action) {
        init(v, op, mapping, composition, identity_monoid, identity_action);
    }
    void init(int n,
              const FuncOperator op,
              const FuncMapping mapping,
              const FuncComposition composition,
              const Monoid &identity_monoid,
              const Action &identity_action) {
        N = n, OP = op, MAPPING = mapping, COMPOSITION = composition;
        IDENTITY_MONOID = identity_monoid, IDENTITY_ACTION = identity_action;
        log = 0, offset = 1;
        while (offset < N) ++log, offset <<= 1;
        dat.assign(offset * 2, IDENTITY_MONOID);
        lazy.assign(offset * 2, IDENTITY_ACTION);
    }
    void init(const vector<Monoid> &v,
              const FuncOperator op,
              const FuncMapping mapping,
              const FuncComposition composition,
              const Monoid &identity_monoid,
              const Action &identity_action) {
        init((int)v.size(), op, mapping, composition, identity_monoid, identity_action);
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
        dat[k] = MAPPING(f, dat[k]);
        if (k < offset) lazy[k] = COMPOSITION(f, lazy[k]);
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
        dat[k] = MAPPING(f, dat[k]);
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
    
    // get max r such that f(v) = True (v = prod(l, r)), O(log N)
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



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// RMQ and RUQ (change the val / min)
void AOJ_RMQ_RUQ() {
    int N, Q;
    cin >> N >> Q;
    const long long identity_monoid = (1LL << 31) - 1;
    const long long identity_action = -1;
    auto op = [&](long long x, long long y) { return min(x, y); };
    auto mapping = [&](long long f, long long x) { return (f != identity_action ? f : x); };
    auto composition = [&](long long g, long long f) { return (g != identity_action ? g : f); };
    LazySegmentTree<long long, long long> seg(N, op, mapping, composition,
                                              identity_monoid, identity_action);
    while (Q--) {
        int type, s, t;
        cin >> type >> s >> t;
        ++t;
        if (type == 0) {
            long long x;
            cin >> x;
            seg.apply(s, t, x);
        } else {
            cout << seg.prod(s, t) << endl;
        }
    }
}

// RSQ and RAQ (add / sum)
void AOJ_RSQ_RAQ() {
    using pll = pair<long long, long long>;
    int N, Q;
    cin >> N >> Q;
    vector<pll> v(N, pll(0, 1));
    const pll identity_monoid = {0, 0};  // {val, range}
    const long long identity_action = 0;
    auto op = [&](pll x, pll y) { return pll(x.first + y.first, x.second + y.second); };
    auto mapping = [&](long long f, pll x) { return pll(x.first + f * x.second, x.second); };
    auto composition = [&](long long g, long long f) { return g + f; };
    LazySegmentTree<pll, long long> seg(v, op, mapping, composition,
                                        identity_monoid, identity_action);
    while (Q--) {
        int type, s, t;
        cin >> type >> s >> t;
        --s;
        if (type == 0) {
            long long x;
            cin >> x;
            seg.apply(s, t, x);
        } else {
            cout << seg.prod(s, t).first << endl;
        }
        //seg.dump();
    }
}

// RMQ and RAQ (add / min)
void AOJ_RMQ_RAQ() {
    int N, Q;
    cin >> N >> Q;
    vector<long long> v(N, 0);
    const long long identity_monoid = (1LL << 31) - 1;
    const long long identity_action = 0;
    auto op = [&](long long x, long long y) { return min(x, y); };
    auto mapping = [&](long long f, long long x) { return f + x; };
    auto composition = [&](long long g, long long f) { return g + f; };
    LazySegmentTree<long long, long long> seg(v, op, mapping, composition,
                                              identity_monoid, identity_action);
    while (Q--) {
        int type, s, t;
        cin >> type >> s >> t;
        ++t;
        if (type == 0) {
            long long x;
            cin >> x;
            seg.apply(s, t, x);
        } else {
            cout << seg.prod(s, t) << endl;
        }
    }
}

// RSQ and RUQ (change the val / sum)
void AOJ_RSQ_RUQ() {
    using pll = pair<long long, long long>;
    int N, Q;
    cin >> N >> Q;
    vector<pll> v(N, pll(0, 1));
    const pll identity_monoid = {0, 0};  // {val, range}
    const long long identity_action = 1 << 30;
    auto op = [&](pll x, pll y) { return pll(x.first + y.first, x.second + y.second); };
    auto mapping = [&](long long f, pll x) {
        return (f != identity_action ? pll(f * x.second, x.second) : x);
    };
    auto composition = [&](long long g, long long f) { return g; };
    LazySegmentTree<pll, long long> seg(v, op, mapping, composition,
                                        identity_monoid, identity_action);
    while (Q--) {
        int type, s, t;
        cin >> type >> s >> t;
        ++t;
        if (type == 0) {
            long long x;
            cin >> x;
            seg.apply(s, t, x);
        } else {
            cout << seg.prod(s, t).first << endl;
        }
    }
}

// modint
template<int MOD> struct Fp {
    long long val;
    constexpr Fp(long long v = 0) noexcept : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr int getmod() const { return MOD; }
    constexpr Fp operator - () const noexcept {
        return val ? MOD - val : 0;
    }
    constexpr Fp operator + (const Fp& r) const noexcept { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp& r) const noexcept { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp& r) const noexcept { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp& r) const noexcept { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp& r) noexcept {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -= (const Fp& r) noexcept {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp& operator *= (const Fp& r) noexcept {
        val = val * r.val % MOD;
        return *this;
    }
    constexpr Fp& operator /= (const Fp& r) noexcept {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        val = val * u % MOD;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr bool operator == (const Fp& r) const noexcept {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp& r) const noexcept {
        return this->val != r.val;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD>& x) noexcept {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD>& x) noexcept {
        return os << x.val;
    }
    friend constexpr Fp<MOD> modpow(const Fp<MOD> &a, long long n) noexcept {
        if (n == 0) return 1;
        auto t = modpow(a, n / 2);
        t = t * t;
        if (n & 1) t = t * a;
        return t;
    }
};

void ACL_Beginner_Contest_E() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    using Node = pair<mint, int>; // val, num
    
    int N, Q;
    cin >> N >> Q;
    vector<mint> ten(N, 1), sum(N+1, 0);
    for (int i = 1; i < N; ++i) ten[i] = ten[i-1] * 10;
    for (int i = 0; i < N; ++i) sum[i+1] = sum[i] + ten[i];

    // define segtree
    auto op = [&](Node x, Node y) {
        mint first = x.first * ten[y.second] + y.first;
        int second = x.second + y.second;
        return Node(first, second);
    };
    auto mapping = [&](int f, Node x) {
        if (f == 0) return x;
        return Node(sum[x.second] * f, x.second);
    };
    auto composition = [&](int g, int f) {
        return g;
    };
    vector<Node> ini(N, Node(mint(1), 1));
    Node identity_monoid = Node(mint(0), 0);
    int identity_action = 0;
    LazySegmentTree<Node, int> seg(ini, op, mapping, composition,
                                   identity_monoid, identity_action);

    // query
    while (Q--) {
        int l, r, d;
        cin >> l >> r >> d;
        --l;
        seg.apply(l, r, d);
        cout << seg.prod(0, N).first << endl;
    }
}


// ABC 322 F - Vacation Query
void ABC_322_F() {
    struct Node {
        int left_zero, left_one, right_zero, right_one, zero, one;
        Node() {}
        Node(int lz, int lo, int rz, int ro, int z, int o)
        : left_zero(lz), left_one(lo), right_zero(rz), right_one(ro), zero(z), one(o) {}
    };
    
    int N, Q;
    string S;
    cin >> N >> Q >> S;
    
    vector<Node> ini(N);
    for (int i = 0; i < N; ++i) {
        if (S[i] == '0') ini[i] = Node(1, 0, 1, 0, 1, 0);
        else ini[i] = Node(0, 1, 0, 1, 0, 1);
    }
    
    Node identity_monoid = Node(-1, -1, -1, -1, -1, -1);
    int identity_action = 0;
    auto op = [&](Node x, Node y) -> Node {
        if (x.one == -1) return y;
        if (y.one == -1) return x;
        Node res;
        res.left_zero = (x.one ? x.left_zero : x.zero + y.left_zero);
        res.left_one = (x.zero ? x.left_one : x.one + y.left_one);
        res.right_zero = (y.one ? y.right_zero : x.right_zero + y.zero);
        res.right_one = (y.zero ? y.right_one : x.right_one + y.one);
        res.zero = max({x.zero, x.right_zero + y.left_zero, y.zero});
        res.one = max({x.one, x.right_one + y.left_one, y.one});
        return res;
    };
    auto mapping = [&](int f, Node x) -> Node {
        if (f) return Node(x.left_one, x.left_zero, x.right_one, x.right_zero, x.one, x.zero);
        else return x;
    };
    auto composition = [&](int g, int f) -> int {
        if (g) return !f;
        else return f;
    };
    
    LazySegmentTree<Node,int> seg(ini, op, mapping, composition,
                                  identity_monoid, identity_action);
    while (Q--) {
        int c, L, R;
        cin >> c >> L >> R;
        --L;
        if (c == 1) {
            seg.apply(L, R, 1);
        } else {
            cout << seg.prod(L, R).one << endl;
        }
    }
}


// ABC 330 E - Mex and Update
void ABC_330_E() {
    int N, Q;
    cin >> N >> Q;
    vector<int> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    
    // 区間加算、区間最小値取得の遅延評価セグメント木
    const int MAX = N + 1;
    vector<int> zero(MAX, 0);
    const int identity_monoid = MAX;
    const int identity_action = 0;
    auto op = [&](int x, int y) { return min(x, y); };
    auto mapping = [&](int f, int x) { return x + f; };
    auto composition = [&](int g, int f) { return g + f; };
    LazySegmentTree<int, int> seg(zero, op, mapping, composition,
                                        identity_monoid, identity_action);
    
    for (int i = 0; i < N; ++i) {
        if (A[i] < MAX) seg.apply(A[i], A[i] + 1, 1);
    }
    
    while (Q--) {
        int id, v;
        cin >> id >> v;
        --id;
        
        // セグ木の更新
        if (A[id] < MAX) seg.apply(A[id], A[id] + 1, -1);
        A[id] = v;
        if (A[id] < MAX) seg.apply(A[id], A[id] + 1, 1);
        
        // セグ木上の二分探索
        // seg.prod(0, r) > 0 を満たす最大の r を求める
        auto check = [&](int val) -> bool { return val > 0; };
        int res = seg.max_right(check, 0);
        cout << res << endl;
    }
}


int main() {
    //AOJ_RMQ_RUQ();
    //AOJ_RSQ_RAQ();
    //AOJ_RMQ_RAQ();
    //AOJ_RSQ_RUQ();
    //ACL_Beginner_Contest_E();
    //ABC_322_F();
    ABC_330_E();
}
