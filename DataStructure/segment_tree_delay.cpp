//
// segment tree (with delay)
//
// verified:
//   AOJ Course Range Query - RMQ and RAQ
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_H&lang=jp
//
//   AtCoder ACL Beginner Contest E - Replace Digits
//     https://atcoder.jp/contests/abl/tasks/abl_e
//


/*
    セグメントツリーは「作用つきモノイド」上で定義される
      FM(a, b): 2 つのモノイド間に定義される演算
      FA(a, d): モノイド元 a への作用素 d による作用
      FC(d, e): 作用素 d への作用素 e の合成
      IDENTITY_MONOID: モノイドの単位元
      IDENTITY_LAZY: 作用素の単位元


    // Construction
    SegTree(N, fm, fa, fc, identity_monoid, identity_lazy)
      ex: starry sky tree (区間加算、区間min取得)
        auto fm = [](long long a, long long b) { return min(a, b); };
        auto fa = [](long long &a, long long d) { a += d; };
        auto fc = [](long long &d, long long e) { d += e; };
        SegTree<long long, long long> seg(N, fm, fa, fc, (1LL<<60), 0);


    // Initialization
    init(n): サイズ n に初期化
    set(a, v): a 番目の値を v にセットする
    build(): set した値を元にセグメントツリー全体を構築する、O(n)


    // Queries
    update(a, b, v): 区間 [a, b) を作用素 v を用いて更新する, O(log n)
    get(a, b): 区間 [a, b) についての演算結果を返す, O(log n)
*/


#include <bits/stdc++.h>
using namespace std;
template<class T1, class T2> ostream& operator << 
(ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }


// Segment Tree
template<class Monoid, class Action> struct SegTree {
    using FuncMonoid = function< Monoid(Monoid, Monoid) >;
    using FuncAction = function< void(Monoid&, Action) >;
    using FuncComposition = function< void(Action&, Action) >;
    FuncMonoid FM;
    FuncAction FA;
    FuncComposition FC;
    Monoid IDENTITY_MONOID;
    Action IDENTITY_LAZY;
    int SIZE, HEIGHT;
    vector<Monoid> dat;
    vector<Action> lazy;
    
    SegTree() {}
    SegTree(int n, const FuncMonoid fm, const FuncAction fa,
            const FuncComposition fc,
            const Monoid &identity_monoid, const Action &identity_lazy)
    : FM(fm), FA(fa), FC(fc), 
      IDENTITY_MONOID(identity_monoid), IDENTITY_LAZY(identity_lazy) {
        SIZE = 1, HEIGHT = 0;
        while (SIZE < n) SIZE <<= 1, ++HEIGHT;
        dat.assign(SIZE * 2, IDENTITY_MONOID);
        lazy.assign(SIZE * 2, IDENTITY_LAZY);
    }
    void init(int n, const FuncMonoid fm, const FuncAction fa,
              const FuncComposition fc,
              const Monoid &identity_monoid, const Action &identity_lazy) {
        FM = fm, FA = fa, FC = fc;
        IDENTITY_MONOID = identity_monoid, IDENTITY_LAZY = identity_lazy;
        SIZE = 1, HEIGHT = 0;
        while (SIZE < n) SIZE <<= 1, ++HEIGHT;
        dat.assign(SIZE * 2, IDENTITY_MONOID);
        lazy.assign(SIZE * 2, IDENTITY_LAZY);
    }
    
    // set, a is 0-indexed
    void set(int a, const Monoid &v) { dat[a + SIZE] = v; }
    void build() {
        for (int k = SIZE - 1; k > 0; --k)
            dat[k] = FM(dat[k*2], dat[k*2+1]);
    }
    
    // update [a, b)
    inline void evaluate(int k) {
        if (lazy[k] == IDENTITY_LAZY) return;
        if (k < SIZE) FC(lazy[k*2], lazy[k]), FC(lazy[k*2+1], lazy[k]);
        FA(dat[k], lazy[k]);
        lazy[k] = IDENTITY_LAZY;
    }
    inline void update(int a, int b, const Action &v, int k, int l, int r) {
        evaluate(k);
        if (a <= l && r <= b) FC(lazy[k], v), evaluate(k);
        else if (a < r && l < b) {
            update(a, b, v, k*2, l, (l+r)>>1);
            update(a, b, v, k*2+1, (l+r)>>1, r);
            dat[k] = FM(dat[k*2], dat[k*2+1]);
        }
    }
    inline void update(int a, int b, const Action &v) { 
        update(a, b, v, 1, 0, SIZE);
    }
    
    // get [a, b)
    inline Monoid get(int a, int b, int k, int l, int r) {
        evaluate(k);
        if (a <= l && r <= b)
            return dat[k];
        else if (a < r && l < b)
            return FM(get(a, b, k*2, l, (l+r)>>1), 
                      get(a, b, k*2+1, (l+r)>>1, r));
        else
            return IDENTITY_MONOID;
    }
    inline Monoid get(int a, int b) { 
        return get(a, b, 1, 0, SIZE);
    }
    inline Monoid operator [] (int a) {
        return get(a, a + 1);
    }
    
    // debug
    void print() {
        for (int i = 0; i < SIZE; ++i) {
            if (i) cout << ",";
            cout << get(i, i+1);
        }
        cout << endl;
    }
};


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

const int MOD = 998244353;
using mint = Fp<MOD>;
using pll = pair<mint,int>; // val, num

int main() {
    int N, Q;
    cin >> N >> Q;
    vector<mint> ten(N, 1), sum(N+1, 0);
    for (int i = 1; i < N; ++i) ten[i] = ten[i-1] * 10;
    for (int i = 0; i < N; ++i) sum[i+1] = sum[i] + ten[i];

    // define segtree
    auto fm = [&](pll a, pll b) {
        mint first = a.first * ten[b.second] + b.first;
        int second = a.second + b.second;
        return pll(first, second);
    };
    auto fa = [&](pll &a, int d) {
        if (d == 0) return;
        a.first = sum[a.second] * d;
    };
    auto fc = [&](int &d, int e) {
        d = e;
    };
    pll identity_monoid = pll(mint(0LL), 0);
    int identity_lazy = 0;
    SegTree<pll, int> seg(N, fm, fa, fc, identity_monoid, identity_lazy);

    // initialization
    for (int i = 0; i < N; ++i) seg.set(i, pll(mint(1LL), 1));
    seg.build();

    // query
    while (Q--) {
        int l, r, d;
        cin >> l >> r >> d;
        --l;
        seg.update(l, r, d);
        cout << seg.get(0, N).first << endl;

        //seg.print();
    }
}
