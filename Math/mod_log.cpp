//
// 離散対数
//   a^x ≡ b (mod m) を満たす非負整数をすべて求める
//　　・唯一解の場合：(x, 0) を返す
//　　・周期的な解の場合：(最小の 0 以上の x, 周期) を返す
//
// verified:
//   Yosupo Library Checker - 	Discrete Logarithm
//     https://judge.yosupo.jp/problem/discrete_logarithm_mod
//
//   AtCoder ABC 222 G - 222
//     https://atcoder.jp/contests/abc222/tasks/abc222_g
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Utility
//------------------------------//

// kth root
// N < 2^64, K <= 64
uint64_t kth_root(uint64_t N, uint64_t K = 2) {
    assert(K >= 1);
    if (N <= 1 || K == 1) return N;
    if (K >= 64) return 1;
    if (N == uint64_t(-1)) --N;
    
    auto mul = [&](uint64_t x, uint64_t y) -> uint64_t {
        if (x < UINT_MAX && y < UINT_MAX) return x * y;
        if (x == uint64_t(-1) || y == uint64_t(-1)) return uint64_t(-1);
        return (x <= uint64_t(-1) / y ? x * y : uint64_t(-1));
    };
    auto power = [&](uint64_t x, uint64_t k) -> uint64_t {
        if (k == 0) return 1ULL;
        uint64_t res = 1ULL;
        while (k) {
            if (k & 1) res = mul(res, x);
            x = mul(x, x);
            k >>= 1;
        }
        return res;
    };
    
    uint64_t res;
    if (K == 2) res = sqrtl(N) - 1;
    else if (K == 3) res = cbrt(N) - 1;
    else res = pow(N, nextafter(1 / double(K), 0));
    while (power(res + 1, K) <= N) ++res;
    return res;
}

// Associative Array
template<class Key, class Val, uint32_t N = 20> struct FastMap {
    static constexpr uint32_t SIZE  = 1u << N;
    static constexpr uint32_t MASK  = SIZE - 1;
    static constexpr uint32_t SHIFT = 64 - N;
    static constexpr uint64_t MAGIC = 11995408973635179863ULL;
    static constexpr uint32_t get_hash(const Key& k) { return k * MAGIC >> SHIFT; }

    // inner values
    array<Key, SIZE> key;
    array<Val, SIZE> val;
    bitset<SIZE> used;

    // constructors
    FastMap() { clear(); }
    void clear() { used.reset(); }
    FastMap(const FastMap&) = default;
    FastMap& operator = (const FastMap&) = default;

    // getters
    Val &operator [] (const Key &k) {
        auto hash = get_hash(k);
        while (true) {
            if (!used[hash]) {
                used[hash] = 1;
                key[hash] = k;
                return val[hash] = Val();
            }
            if (key[hash] == k) return val[hash];
            ++hash &= MASK;
        }
    }
    const Val &operator [] (const Key &k) const {
        auto hash = get_hash(k);
        while (true) {
            if (!used[hash]) {
                used[hash] = 1;
                key[hash] = k;
                return val[hash] = Val();
            }
            if (key[hash] == k) return val[hash];
            ++hash &= MASK;
        }
    }
    Val get(const Key &k) const {
        auto hash = get_hash(k);
        while (true) {
            if (!used[hash]) return Val();
            if (key[hash] == k) return val[hash];
            ++hash &= MASK;
        }
    }
    bool count(const Key &k) const {
        auto hash = get_hash(k);
        while (true) {
            if (!used[hash]) return false;
            if (key[hash] == k) return true;
            ++hash &= MASK;
        }
    }
};

//-----------------------------------------//
// necessary number-theoretic algorithms
//-----------------------------------------//

// mod pow
template<class T_VAL, class T_MOD>
constexpr T_VAL mod_pow(T_VAL a, T_VAL n, T_MOD m) {
    T_VAL res = 1;
    while (n > 0) {
        if (n % 2 == 1) res = res * a % m;
        a = a * a % m;
        n >>= 1;
    }
    return res;
}

// montgomery modint (MOD < 2^62, MOD is odd)
struct MontgomeryModInt64 {
    using mint = MontgomeryModInt64;
    using u64 = uint64_t;
    using u128 = __uint128_t;
    
    // static menber
    static u64 MOD;
    static u64 INV_MOD;  // INV_MOD * MOD ≡ 1 (mod 2^64)
    static u64 T128;  // 2^128 (mod MOD)
    
    // inner value
    u64 val;
    
    // constructor
    MontgomeryModInt64() : val(0) { }
    MontgomeryModInt64(long long v) : val(reduce((u128(v) + MOD) * T128)) { }
    u64 get() const {
        u64 res = reduce(val);
        return res >= MOD ? res - MOD : res;
    }
    
    // mod getter and setter
    static u64 get_mod() { return MOD; }
    static void set_mod(u64 mod) {
        assert(mod < (1LL << 62));
        assert((mod & 1));
        MOD = mod;
        T128 = -u128(mod) % mod;
        INV_MOD = get_inv_mod();
    }
    static u64 get_inv_mod() {
        u64 res = MOD;
        for (int i = 0; i < 5; ++i) res *= 2 - MOD * res;
        return res;
    }
    static u64 reduce(const u128 &v) {
        return (v + u128(u64(v) * u64(-INV_MOD)) * MOD) >> 64;
    }
    
    // arithmetic operators
    mint operator + () const { return mint(*this); }
    mint operator - () const { return mint() - mint(*this); }
    mint operator + (const mint &r) const { return mint(*this) += r; }
    mint operator - (const mint &r) const { return mint(*this) -= r; }
    mint operator * (const mint &r) const { return mint(*this) *= r; }
    mint operator / (const mint &r) const { return mint(*this) /= r; }
    mint& operator += (const mint &r) {
        if ((val += r.val) >= 2 * MOD) val -= 2 * MOD;
        return *this;
    }
    mint& operator -= (const mint &r) {
        if ((val += 2 * MOD - r.val) >= 2 * MOD) val -= 2 * MOD;
        return *this;
    }
    mint& operator *= (const mint &r) {
        val = reduce(u128(val) * r.val);
        return *this;
    }
    mint& operator /= (const mint &r) {
        *this *= r.inv();
        return *this;
    }
    mint inv() const { return pow(MOD - 2); }
    mint pow(u128 n) const {
        mint res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }

    // other operators
    bool operator == (const mint &r) const {
        return (val >= MOD ? val - MOD : val) == (r.val >= MOD ? r.val - MOD : r.val);
    }
    bool operator != (const mint &r) const {
        return (val >= MOD ? val - MOD : val) != (r.val >= MOD ? r.val - MOD : r.val);
    }
    mint& operator ++ () {
        ++val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    mint& operator -- () {
        if (val == 0) val += MOD;
        --val;
        return *this;
    }
    mint operator ++ (int) {
        mint res = *this;
        ++*this;
        return res;
    }
    mint operator -- (int) {
        mint res = *this;
        --*this;
        return res;
    }
    friend istream& operator >> (istream &is, mint &x) {
        long long t;
        is >> t;
        x = mint(t);
        return is;
    }
    friend ostream& operator << (ostream &os, const mint &x) {
        return os << x.get();
    }
    friend mint pow(const mint &r, long long n) {
        return r.pow(n);
    }
    friend mint inv(const mint &r) {
        return r.inv();
    }
};

typename MontgomeryModInt64::u64
MontgomeryModInt64::MOD, MontgomeryModInt64::INV_MOD, MontgomeryModInt64::T128;

// Miller-Rabin
bool MillerRabin(long long N, const vector<long long> &A) {
    assert(N % 2 == 1);
    assert(N < (1LL<<62));
    using mint = MontgomeryModInt64;
    mint::set_mod(N);
    
    long long s = 0, d = N - 1;
    while (d % 2 == 0) {
        ++s;
        d >>= 1;
    }
    for (auto a : A) {
        if (N <= a) return true;
        mint x = mint(a).pow(d);
        if (x != 1) {
            long long t;
            for (t = 0; t < s; ++t) {
                if (x == N - 1) break;
                x *= x;
            }
            if (t == s) return false;
        }
    }
    return true;
}

bool is_prime(long long N) {
    if (N <= 1) return false;
    else if (N == 2) return true;
    else if (N % 2 == 0) return false;
    else if (N < 4759123141LL)
        return MillerRabin(N, {2, 7, 61});
    else
        return MillerRabin(N, {2, 325, 9375, 28178, 450775, 9780504, 1795265022});
}

// Pollard's Rho
unsigned int xor_shift_rng() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty, ty = tz, tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}

long long pollard(long long N) {
    if (N % 2 == 0) return 2;
    if (is_prime(N)) return N;
    
    assert(N < (1LL<<62));
    using mint = MontgomeryModInt64;
    mint::set_mod(N);
    
    long long step = 0;
    while (true) {
        mint r = xor_shift_rng();  // random r
        auto f = [&](mint x) -> mint { return x * x + r; };
        mint x = ++step, y = f(x);
        while (true) {
            long long p = gcd((y - x).get(), N);
            if (p == 0 || p == N) break;
            if (p != 1) return p;
            x = f(x);
            y = f(f(y));
        }
    }
}

vector<long long> pollard_prime_factorize(long long N) {
    if (N == 1) return {};
    long long p = pollard(N);
    if (p == N) return {p};
    vector<long long> left = pollard_prime_factorize(p);
    vector<long long> right = pollard_prime_factorize(N / p);
    if (left.size() > right.size()) swap(left, right);
    left.insert(left.end(), right.begin(), right.end());
    sort(left.begin(), left.end());
    return left;
}

vector<pair<long long, long long>> prime_factorize(long long N) {
    vector<pair<long long, long long>> res;
    const auto &prs = pollard_prime_factorize(N);
    long long prev = -1, num = 0;
    for (const auto &pr : prs) {
        if (pr == prev) ++num;
        else {
            if (prev != -1) res.emplace_back(prev, num);
            prev = pr, num = 1;
        }
    }
    if (prev != -1) res.emplace_back(prev, num);
    return res;
}

// ax + by = gcd(a, b)
// return: gcd(a, b)
template<class T_VAL> T_VAL ext_gcd(T_VAL a, T_VAL b, T_VAL &x, T_VAL &y) {
    T_VAL xsign = 1, ysign = 1;
    if (a < 0) {
        a = -a, xsign *= -1;
    }
    if (b < 0) {
        b = -b, ysign *= -1;
    }
    
    if (b == 0) {
        x = xsign, y = 0;
        return a;
    }
    T_VAL d = ext_gcd(b, a % b, y, x);
    y -= a / b * x;
    x *= xsign, y *= ysign;
    return d;
}


//------------------------------//
// discrete logarithm
//------------------------------//

// calc order
template<class T_VAL, class T_MOD>
constexpr T_VAL calc_order(T_VAL a, T_MOD m) {
    assert(gcd(a, m) == 1);
    T_VAL euler = m;
    for (const auto &[p, exp] : prime_factorize(m)) euler -= euler / p;
    T_VAL res = euler;
    for (const auto &[p, exp] : prime_factorize(euler)) {
        while (res % p == 0 && mod_pow(a, res / p, m) == 1) res /= p;
    }
    return res;
}

// Discrete Logarithm
// find min positive x s.t. a^x = b (mod m)
// return (c, d) where solution is x = c (mod d)
template<class T_VAL, class T_MOD>
constexpr pair<T_VAL, T_MOD> mod_log(T_VAL a, T_VAL b, T_MOD m) {
    if (m == 1) return {0, 1};
    if (a == 0 && b == 0) return {1, 1};
    T_VAL m1 = m, pw = 1, div = gcd(a, m1), res = 0, x, y;
    for (; (div = gcd(a, m1)) > 1; res++, m1 /= div, pw = pw * a % m) {
        if (pw == b) return {res, 0};  // aperiodic solution
    }
    auto g = ext_gcd(pw, m, x, y);
    if (b % g > 0) return {-1, 0};  // no solution
    b = x * (b / g) % m;
    if (b < 0) b += m;
    T_VAL order = calc_order(a, m1), sq = kth_root(order, 2) + 1, an = mod_pow(a, sq, m1);
    FastMap<T_VAL, T_VAL> half;
    pw = an;
    for (T_VAL p = 1; p <= sq; p++, pw = pw * an % m1) half[pw] = p;
    pw = b % m1;
    for (T_VAL q = 0; q <= sq; q++, pw = pw * a % m1) {
        if (half.count(pw)) {
            res += (sq * half[pw] - q) % order;
            return {res, order};  // periodic solution
        }
    }
    return {-1, 0};  // no solution
}


//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Discrete Logarithm
void Yosupo_Discrete_Logarithm() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int T;
    cin >> T;
    while (T--) {
        long long x, y, m;
        cin >> x >> y >> m;
        auto [exp, syuki] = mod_log(x, y, m);
        cout << exp << '\n';
    }
}

// AtCoder ABC 222 G - 222
void ABC_222_G() {
    int T;
    cin >> T;
    while (T--) {
        long long K;
        cin >> K;
        if (K % 2 == 0) K /= 2, K *= 9;
        else K *= 9;
        if (K % 2 == 0 || K % 5 == 0) cout << -1 << endl;
        else {
            auto [exp, syuki] = mod_log(10LL, 1LL, K);
            cout << syuki << endl;
        }
    }
}


int main() {
    Yosupo_Discrete_Logarithm();
    //ABC_222_G();
}