//
// SWAG (Sliding Window Aggregation)
//
// references:
//   https://scrapbox.io/data-structures/Sliding_Window_Aggregation
//   https://www.slideshare.net/catupper/amortize-analysis-of-deque-with-2-stack
//
// verified:
//   AOJ DSL_3_D - Sliding Minimum Elements
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_3_D&lang=ja
//
//   Yosupo Library Checker - Deque Operate All Composite
//     https://judge.yosupo.jp/problem/deque_operate_all_composite
//


#include <bits/stdc++.h>
using namespace std;


// SWAG
template<class Monoid> struct SWAG {
    using Func = function<Monoid(Monoid, Monoid)>;
    
    // core member
    Func OP;
    Monoid IDENTITY;
    
    // inner data
    int siz;
    vector<Monoid> dat_left, dat_right, sum_left, sum_right;
    
    // constructor
    SWAG() {}
    SWAG(const Func &op, const Monoid &identity) {
        init(op, identity);
    }
    SWAG(const vector<Monoid> &vec, const Func &op, const Monoid &identity) {
        init(vec, op, identity);
    }
    void init(const Func &op, const Monoid &identity) {
        OP = op;
        IDENTITY = identity;
        clear();
    }
    void init(const vector<Monoid> &vec, const Func &op, const Monoid &identity) {
        init(op, identity);
        for (const auto &v : vec) push_back(v);
    }
    void clear() {
        siz = 0;
        dat_left.clear(), dat_right.clear();
        sum_left = {IDENTITY}, sum_right = {IDENTITY};
    }
    
    // getter
    int size() { return siz; }
    
    // push
    void push_back(const Monoid &v) {
        ++siz;
        dat_right.emplace_back(v);
        sum_right.emplace_back(OP(sum_right.back(), v));
    }
    void push_front(const Monoid &v) {
        ++siz;
        dat_left.emplace_back(v);
        sum_left.emplace_back(OP(v, sum_left.back()));
    }
    
    // pop
    void rebuild() {
        vector<Monoid> tmp;
        for (int i = dat_left.size() - 1; i >= 0; --i) tmp.emplace_back(dat_left[i]);
        for (int i = 0; i < dat_right.size(); ++i) tmp.emplace_back(dat_right[i]);
        //tmp.insert(tmp.end(), dat_right.begin(), dat_right.end());
        clear();
        int mid = tmp.size() / 2;
        for (int i = mid - 1; i >= 0; --i) push_front(tmp[i]);
        for (int i = mid; i < tmp.size(); ++i) push_back(tmp[i]);
        assert(siz == tmp.size());
    }
    void pop_back() {
        if (siz == 1) return clear();
        if (dat_right.empty()) rebuild();
        --siz;
        dat_right.pop_back();
        sum_right.pop_back();
    }
    void pop_front() {
        if (siz == 1) return clear();
        if (dat_left.empty()) rebuild();
        --siz;
        dat_left.pop_back();
        sum_left.pop_back();
    }
    
    // prod
    Monoid prod() {
        return OP(sum_left.back(), sum_right.back());
    }
    
    // debug
    friend ostream& operator << (ostream &s, const SWAG &sw) {
        for (int i = sw.dat_left.size() - 1; i >= 0; --i) s << sw.dat_left[i] << " ";
        for (int i = 0; i < sw.dat_right.size(); ++i) s << sw.dat_right[i] << " ";
        return s;
    }
};



//------------------------------//
// Examples
//------------------------------//

// AOJ DSL_3_D - Sliding Minimum Elements
void AOJ_DSL_3_D() {
    int N, L;
    cin >> N >> L;
    vector<int> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    
    const int INF = 1<<30;
    SWAG<int> sw([&](int l, int r){return min(l, r);}, INF);
    for (int i = 0; i < L; ++i) sw.push_back(a[i]);

    for (int i = L; i < N; ++i) {
        cout << sw.prod() << " ";
        sw.pop_front();
        sw.push_back(a[i]);
    }
    cout << sw.prod() << endl;
}


// Yosupo Library Checker - Deque Operate All Composite
// modint
template<int MOD> struct Fp {
    // inner value
    long long val;
    
    // constructor
    constexpr Fp() : val(0) { }
    constexpr Fp(long long v) : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr long long get() const { return val; }
    constexpr int get_mod() const { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator + () const { return Fp(*this); }
    constexpr Fp operator - () const { return Fp(0) - Fp(*this); }
    constexpr Fp operator + (const Fp &r) const { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp &r) {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) {
        val = val * r.val % MOD;
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) {
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
    constexpr Fp pow(long long n) const {
        Fp res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const {
        Fp res(1), div(*this);
        return res / div;
    }

    // other operators
    constexpr bool operator == (const Fp &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const {
        return this->val != r.val;
    }
    constexpr Fp& operator ++ () {
        ++val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -- () {
        if (val == 0) val += MOD;
        --val;
        return *this;
    }
    constexpr Fp operator ++ (int) const {
        Fp res = *this;
        ++*this;
        return res;
    }
    constexpr Fp operator -- (int) const {
        Fp res = *this;
        --*this;
        return res;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD> &x) {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD> &x) {
        return os << x.val;
    }
    friend constexpr Fp<MOD> pow(const Fp<MOD> &r, long long n) {
        return r.pow(n);
    }
    friend constexpr Fp<MOD> inv(const Fp<MOD> &r) {
        return r.inv();
    }
};
void yosupo_deque_operate_all_composite() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    using Monoid = pair<mint,mint>;
    auto op = [&](Monoid x, Monoid y) {
        return Monoid(x.first * y.first, x.second * y.first + y.second);
    };
    Monoid identity = {1, 0};
    
    SWAG<Monoid> sw(op, identity);
    int Q;
    cin >> Q;
    while (Q--) {
        int t, a, b, x;
        cin >> t;
        if (t == 0) {
            cin >> a >> b;
            sw.push_front(Monoid(a, b));
        } else if (t == 1) {
            cin >> a >> b;
            sw.push_back(Monoid(a, b));
        } else if (t == 2) {
            sw.pop_front();
        } else if (t == 3) {
            sw.pop_back();
        } else {
            cin >> x;
            auto f = sw.prod();
            cout << f.first * x + f.second << endl;
        }
    }
}


int main() {
    //AOJ_DSL_3_D();
    yosupo_deque_operate_all_composite();
}

