//
// Formal Power Series
//
// verified:
//   Yosupo Library Checker - Inv of Formal Power Series
//     https://judge.yosupo.jp/problem/inv_of_formal_power_series
//
//   Yosupo Library Checker - Exp of Formal Power Series
//     https://judge.yosupo.jp/problem/exp_of_formal_power_series
//
//   Yosupo Library Checker - Log of Formal Power Series
//     https://judge.yosupo.jp/problem/log_of_formal_power_series
//
//   Yosupo Library Checker - Pow of Formal Power Series
//     https://judge.yosupo.jp/problem/pow_of_formal_power_series
//
//   Yosupo Library Checker - Sqrt of Formal Power Series
//     https://judge.yosupo.jp/problem/sqrt_of_formal_power_series
//
//   HackerRank Array Restoring
//     https://www.hackerrank.com/contests/happy-query-contest/challenges/array-restoring/problem
//
//   Codeforces 205 Div1 E. The Child and Binary TreeE
//     https://codeforces.com/contest/438/problem/E
//
//   TDPC T - フィボナッチ (mod. 1000000007)
//     https://atcoder.jp/contests/tdpc/tasks/tdpc_fibonacci
//


#include <bits/stdc++.h>
using namespace std;


// modint
template<int MOD> struct Fp {
    // inner value
    long long val;
    
    // constructor
    constexpr Fp() noexcept : val(0) { }
    constexpr Fp(long long v) noexcept : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr long long get() const noexcept { return val; }
    constexpr int get_mod() const noexcept { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator - () const noexcept {
        return val ? MOD - val : 0;
    }
    constexpr Fp operator + (const Fp &r) const noexcept { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const noexcept { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const noexcept { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const noexcept { return Fp(*this) /= r; }
    constexpr Fp& operator ++ () noexcept {
        ++val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -- () noexcept {
        if (val == 0) val += MOD;
        --val;
    }
    constexpr Fp& operator += (const Fp &r) noexcept {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) noexcept {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) noexcept {
        val = val * r.val % MOD;
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) noexcept {
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
    constexpr Fp pow(long long n) const noexcept {
        Fp res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const noexcept {
        Fp res(1), div(*this);
        return res / div;
    }

    // other operators
    constexpr bool operator == (const Fp &r) const noexcept {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const noexcept {
        return this->val != r.val;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD> &x) noexcept {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD> &x) noexcept {
        return os << x.val;
    }
    friend constexpr Fp<MOD> modpow(const Fp<MOD> &r, long long n) noexcept {
        return r.pow(n);
    }
    friend constexpr Fp<MOD> modinv(const Fp<MOD> &r) noexcept {
        return r.inv();
    }
};

namespace NTT {
    long long modpow(long long a, long long n, int mod) {
        long long res = 1;
        while (n > 0) {
            if (n & 1) res = res * a % mod;
            a = a * a % mod;
            n >>= 1;
        }
        return res;
    }

    long long modinv(long long a, int mod) {
        long long b = mod, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        u %= mod;
        if (u < 0) u += mod;
        return u;
    }

    int calc_primitive_root(int mod) {
        if (mod == 2) return 1;
        if (mod == 167772161) return 3;
        if (mod == 469762049) return 3;
        if (mod == 754974721) return 11;
        if (mod == 998244353) return 3;
        int divs[20] = {};
        divs[0] = 2;
        int cnt = 1;
        long long x = (mod - 1) / 2;
        while (x % 2 == 0) x /= 2;
        for (long long i = 3; i * i <= x; i += 2) {
            if (x % i == 0) {
                divs[cnt++] = i;
                while (x % i == 0) x /= i;
            }
        }
        if (x > 1) divs[cnt++] = x;
        for (int g = 2;; g++) {
            bool ok = true;
            for (int i = 0; i < cnt; i++) {
                if (modpow(g, (mod - 1) / divs[i], mod) == 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) return g;
        }
    }

    int get_fft_size(int N, int M) {
        int size_a = 1, size_b = 1;
        while (size_a < N) size_a <<= 1;
        while (size_b < M) size_b <<= 1;
        return max(size_a, size_b) << 1;
    }

    // number-theoretic transform
    template<class mint> void trans(vector<mint> &v, bool inv = false) {
        if (v.empty()) return;
        int N = (int)v.size();
        int MOD = v[0].get_mod();
        int PR = calc_primitive_root(MOD);
        static bool first = true;
        static vector<long long> vbw(30), vibw(30);
        if (first) {
            first = false;
            for (int k = 0; k < 30; ++k) {
                vbw[k] = modpow(PR, (MOD - 1) >> (k + 1), MOD);
                vibw[k] = modinv(vbw[k], MOD);
            }
        }
        for (int i = 0, j = 1; j < N - 1; j++) {
            for (int k = N >> 1; k > (i ^= k); k >>= 1);
            if (i > j) swap(v[i], v[j]);
        }
        for (int k = 0, t = 2; t <= N; ++k, t <<= 1) {
            long long bw = vbw[k];
            if (inv) bw = vibw[k];
            for (int i = 0; i < N; i += t) {
                mint w = 1;
                for (int j = 0; j < t/2; ++j) {
                    int j1 = i + j, j2 = i + j + t/2;
                    mint c1 = v[j1], c2 = v[j2] * w;
                    v[j1] = c1 + c2;
                    v[j2] = c1 - c2;
                    w *= bw;
                }
            }
        }
        if (inv) {
            long long invN = modinv(N, MOD);
            for (int i = 0; i < N; ++i) v[i] = v[i] * invN;
        }
    }

    // for garner
    static constexpr int MOD0 = 754974721;
    static constexpr int MOD1 = 167772161;
    static constexpr int MOD2 = 469762049;
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    static const mint1 imod0 = 95869806; // modinv(MOD0, MOD1);
    static const mint2 imod1 = 104391568; // modinv(MOD1, MOD2);
    static const mint2 imod01 = 187290749; // imod1 / MOD0;

    // small case (T = mint, long long)
    template<class T> vector<T> naive_mul(const vector<T> &A, const vector<T> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        vector<T> res(N + M - 1);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < M; ++j)
                res[i + j] += A[i] * B[j];
        return res;
    }

    // mint
    template<class mint> vector<mint> mul(const vector<mint> &A, const vector<mint> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        if (min(N, M) < 30) return naive_mul(A, B);
        int MOD = A[0].get_mod();
        int size_fft = get_fft_size(N, M);
        if (MOD == 998244353) {
            vector<mint> a(size_fft), b(size_fft), c(size_fft);
            for (int i = 0; i < N; ++i) a[i] = A[i];
            for (int i = 0; i < M; ++i) b[i] = B[i];
            trans(a), trans(b);
            vector<mint> res(size_fft);
            for (int i = 0; i < size_fft; ++i) res[i] = a[i] * b[i];
            trans(res, true);
            res.resize(N + M - 1);
            return res;
        }
        vector<mint0> a0(size_fft, 0), b0(size_fft, 0), c0(size_fft, 0);
        vector<mint1> a1(size_fft, 0), b1(size_fft, 0), c1(size_fft, 0);
        vector<mint2> a2(size_fft, 0), b2(size_fft, 0), c2(size_fft, 0);
        for (int i = 0; i < N; ++i)
            a0[i] = A[i].val, a1[i] = A[i].val, a2[i] = A[i].val;
        for (int i = 0; i < M; ++i)
            b0[i] = B[i].val, b1[i] = B[i].val, b2[i] = B[i].val;
        trans(a0), trans(a1), trans(a2), trans(b0), trans(b1), trans(b2);
        for (int i = 0; i < size_fft; ++i) {
            c0[i] = a0[i] * b0[i];
            c1[i] = a1[i] * b1[i];
            c2[i] = a2[i] * b2[i];
        }
        trans(c0, true), trans(c1, true), trans(c2, true);
        static const mint mod0 = MOD0, mod01 = mod0 * MOD1;
        vector<mint> res(N + M - 1);
        for (int i = 0; i < N + M - 1; ++i) {
            int y0 = c0[i].val;
            int y1 = (imod0 * (c1[i] - y0)).val;
            int y2 = (imod01 * (c2[i] - y0) - imod1 * y1).val;
            res[i] = mod01 * y2 + mod0 * y1 + y0;
        }
        return res;
    }

    // long long
    vector<long long> mul_ll(const vector<long long> &A, const vector<long long> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        if (min(N, M) < 30) return naive_mul(A, B);
        int size_fft = get_fft_size(N, M);
        vector<mint0> a0(size_fft, 0), b0(size_fft, 0), c0(size_fft, 0);
        vector<mint1> a1(size_fft, 0), b1(size_fft, 0), c1(size_fft, 0);
        vector<mint2> a2(size_fft, 0), b2(size_fft, 0), c2(size_fft, 0);
        for (int i = 0; i < N; ++i)
            a0[i] = A[i], a1[i] = A[i], a2[i] = A[i];
        for (int i = 0; i < M; ++i)
            b0[i] = B[i], b1[i] = B[i], b2[i] = B[i];
        trans(a0), trans(a1), trans(a2), trans(b0), trans(b1), trans(b2);
        for (int i = 0; i < size_fft; ++i) {
            c0[i] = a0[i] * b0[i];
            c1[i] = a1[i] * b1[i];
            c2[i] = a2[i] * b2[i];
        }
        trans(c0, true), trans(c1, true), trans(c2, true);
        static const long long mod0 = MOD0, mod01 = mod0 * MOD1;
        vector<long long> res(N + M - 1);
        for (int i = 0; i < N + M - 1; ++i) {
            int y0 = c0[i].val;
            int y1 = (imod0 * (c1[i] - y0)).val;
            int y2 = (imod01 * (c2[i] - y0) - imod1 * y1).val;
            res[i] = mod01 * y2 + mod0 * y1 + y0;
        }
        return res;
    }
};

// Formal Power Series
template<typename mint> struct FPS : vector<mint> {
    using vector<mint>::vector;
 
    // constructor
    constexpr FPS(const vector<mint> &r) : vector<mint>(r) {}
 
    // core operator
    constexpr FPS pre(int siz) const {
        return FPS(begin(*this), begin(*this) + min((int)this->size(), siz));
    }
    constexpr FPS rev() const {
        FPS res = *this;
        reverse(begin(res), end(res));
        return res;
    }
    constexpr FPS& normalize() {
        while (!this->empty() && this->back() == 0) this->pop_back();
        return *this;
    }
 
    // basic operator
    constexpr FPS operator - () const noexcept {
        FPS res = (*this);
        for (int i = 0; i < (int)res.size(); ++i) res[i] = -res[i];
        return res;
    }
    constexpr FPS operator + (const mint &v) const { return FPS(*this) += v; }
    constexpr FPS operator + (const FPS &r) const { return FPS(*this) += r; }
    constexpr FPS operator - (const mint &v) const { return FPS(*this) -= v; }
    constexpr FPS operator - (const FPS &r) const { return FPS(*this) -= r; }
    constexpr FPS operator * (const mint &v) const { return FPS(*this) *= v; }
    constexpr FPS operator * (const FPS &r) const { return FPS(*this) *= r; }
    constexpr FPS operator / (const mint &v) const { return FPS(*this) /= v; }
    constexpr FPS operator / (const FPS &r) const { return FPS(*this) /= r; }
    constexpr FPS operator % (const FPS &r) const { return FPS(*this) %= r; }
    constexpr FPS operator << (int x) const { return FPS(*this) <<= x; }
    constexpr FPS operator >> (int x) const { return FPS(*this) >>= x; }
    constexpr FPS& operator += (const mint &v) {
        if (this->empty()) this->resize(1);
        (*this)[0] += v;
        return *this;
    }
    constexpr FPS& operator += (const FPS &r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return this->normalize();
    }
    constexpr FPS& operator -= (const mint &v) {
        if (this->empty()) this->resize(1);
        (*this)[0] -= v;
        return *this;
    }
    constexpr FPS& operator -= (const FPS &r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return this->normalize();
    }
    constexpr FPS& operator *= (const mint &v) {
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= v;
        return *this;
    }
    constexpr FPS& operator *= (const FPS &r) {
        return *this = NTT::mul((*this), r);
    }
    constexpr FPS& operator /= (const mint &v) {
        assert(v != 0);
        mint iv = modinv(v);
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= iv;
        return *this;
    }
    
    // division, r must be normalized (r.back() must not be 0)
    constexpr FPS& operator /= (const FPS &r) {
        assert(!r.empty());
        assert(r.back() != 0);
        this->normalize();
        if (this->size() < r.size()) {
            this->clear();
            return *this;
        }
        int need = (int)this->size() - (int)r.size() + 1;
        *this = (rev().pre(need) * r.rev().inv(need)).pre(need).rev();
        return *this;
    }
    constexpr FPS& operator %= (const FPS &r) {
        assert(!r.empty());
        assert(r.back() != 0);
        this->normalize();
        FPS q = (*this) / r;
        return *this -= q * r;
    }
    constexpr FPS& operator <<= (int x) {
        FPS res(x, 0);
        res.insert(res.end(), begin(*this), end(*this));
        return *this = res;
    }
    constexpr FPS& operator >>= (int x) {
        FPS res;
        res.insert(res.end(), begin(*this) + x, end(*this));
        return *this = res;
    }
    constexpr mint eval(const mint &v) {
        mint res = 0;
        for (int i = (int)this->size()-1; i >= 0; --i) {
            res *= v;
            res += (*this)[i];
        }
        return res;
    }

    // advanced operation
    // df/dx
    constexpr FPS diff() const {
        int n = (int)this->size();
        FPS res(n-1);
        for (int i = 1; i < n; ++i) res[i-1] = (*this)[i] * i;
        return res;
    }
    
    // \int f dx
    constexpr FPS integral() const {
        int n = (int)this->size();
        FPS res(n+1, 0);
        for (int i = 0; i < n; ++i) res[i+1] = (*this)[i] / (i+1);
        return res;
    }
    
    // inv(f), f[0] must not be 0
    constexpr FPS inv(int deg) const {
        assert((*this)[0] != 0);
        if (deg < 0) deg = (int)this->size();
        FPS res({mint(1) / (*this)[0]});
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + res - res * res * pre(i << 1)).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS inv() const {
        return inv((int)this->size());
    }
    
    // log(f) = \int f'/f dx, f[0] must be 1
    constexpr FPS log(int deg) const {
        assert((*this)[0] == 1);
        FPS res = (diff() * inv(deg)).integral();
        res.resize(deg);
        return res;
    }
    constexpr FPS log() const {
        return log((int)this->size());
    }
    
    // exp(f), f[0] must be 0
    constexpr FPS exp(int deg) const {
        assert((*this)[0] == 0);
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = res * (pre(i << 1) - res.log(i << 1) + 1).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS exp() const {
        return exp((int)this->size());
    }
    
    // pow(f) = exp(e * log f)
    constexpr FPS pow(long long e, int deg) const {
        if (e == 0) {
            FPS res(deg, 0);
            res[0] = 1;
            return res;
        }
        long long i = 0;
        while (i < (int)this->size() && (*this)[i] == 0) ++i;
        if (i == (int)this->size() || i > (deg - 1) / e) return FPS(deg, 0);
        mint k = (*this)[i];
        FPS res = ((((*this) >> i) / k).log(deg) * e).exp(deg) * mint(k).pow(e) << (e * i);
        res.resize(deg);
        return res;
    }
    constexpr FPS pow(long long e) const {
        return pow(e, (int)this->size());
    }
    
    // sqrt(f), f[0] must be 1
    constexpr FPS sqrt_base(int deg) const {
        assert((*this)[0] == 1);
        mint inv2 = mint(1) / 2;
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + pre(i << 1) * res.inv(i << 1)).pre(i << 1);
            for (mint &x : res) x *= inv2;
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS sqrt_base() const {
        return sqrt_base((int)this->size());
    }
    
    // friend operators
    friend constexpr FPS diff(const FPS &f) { return f.diff(); }
    friend constexpr FPS integral(const FPS &f) { return f.integral(); }
    friend constexpr FPS inv(const FPS &f, int deg) { return f.inv(deg); }
    friend constexpr FPS inv(const FPS &f) { return f.inv((int)f.size()); }
    friend constexpr FPS log(const FPS &f, int deg) { return f.log(deg); }
    friend constexpr FPS log(const FPS &f) { return f.log((int)f.size()); }
    friend constexpr FPS exp(const FPS &f, int deg) { return f.exp(deg); }
    friend constexpr FPS exp(const FPS &f) { return f.exp((int)f.size()); }
    friend constexpr FPS pow(const FPS &f, long long e, int deg) { return f.pow(e, deg); }
    friend constexpr FPS pow(const FPS &f, long long e) { return f.pow(e, (int)f.size()); }
    friend constexpr FPS sqrt_base(const FPS &f, int deg) { return f.sqrt_base(deg); }
    friend constexpr FPS sqrt_base(const FPS &f) { return f.sqrt_base((int)f.size()); }
};


/*/////////////////////////////*/
// Polynomial, FPS algorithms
/*/////////////////////////////*/

// Bostan-Mori
// find [x^N] P(x)/Q(x), O(K log K log N)
// deg(Q(x)) = K, deg(P(x)) < K
template<typename mint> mint BostanMori(const FPS<mint> &P, const FPS<mint> &Q, long long N) {
    assert(!P.empty() && !Q.empty());
    if (N == 0 || Q.size() == 1) return P[0] / Q[0];
    
    int qdeg = (int)Q.size();
    FPS<mint> P2{P}, minusQ{Q};
    P2.resize(qdeg - 1);
    for (int i = 1; i < (int)Q.size(); i += 2) minusQ[i] = -minusQ[i];
    P2 *= minusQ;
    FPS<mint> Q2 = Q * minusQ;
    FPS<mint> S(qdeg - 1), T(qdeg);
    for (int i = 0; i < (int)S.size(); ++i) {
        S[i] = (N % 2 == 0 ? P2[i * 2] : P2[i * 2 + 1]);
    }
    for (int i = 0; i < (int)T.size(); ++i) {
        T[i] = Q2[i * 2];
    }
    return BostanMori(S, T, N >> 1);
}



/*/////////////////////////////*/
// solvers
/*/////////////////////////////*/

void Yosupo_Inv_of_FPS() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    int N;
    cin >> N;
    FPS<mint> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];

    auto res = inv(a);
    for (int i = 0; i < res.size(); ++i) {
        if (i) cout << " ";
        cout << res[i];
    }
    cout << endl;
}

void Yosupo_Exp_of_FPS() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    int N;
    cin >> N;
    FPS<mint> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];

    auto res = exp(a);
    for (int i = 0; i < res.size(); ++i) {
        if (i) cout << " ";
        cout << res[i];
    }
    cout << endl;
}

void Yosupo_Log_of_FPS() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    int N;
    cin >> N;
    FPS<mint> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];

    auto res = log(a);
    for (int i = 0; i < res.size(); ++i) {
        if (i) cout << " ";
        cout << res[i];
    }
    cout << endl;
}

void Yosupo_Pow_of_FPS() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    long long N, M;
    cin >> N >> M;
    FPS<mint> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];

    auto res = pow(a, M);
    for (int i = 0; i < res.size(); ++i) {
        if (i) cout << " ";
        cout << res[i];
    }
    cout << endl;
}

void Yosupo_Sqrt_of_FPS() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    int N;
    cin >> N;
    FPS<mint> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];

    auto res = sqrt_base(a);
    for (int i = 0; i < res.size(); ++i) {
        if (i) cout << " ";
        cout << res[i];
    }
    cout << endl;
}

void HackerRankArrayRestoring() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    int N, M, Q;
    cin >> N >> M >> Q;
    FPS<mint> f(N-M+1, 0);
    for (int i = 0; i < Q; ++i) {
        int l, r;
        cin >> l >> r;
        --l;
        f[l] += 1;
    }
    f.normalize();

    FPS<mint> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    auto B = A / f;
    B.resize(M);
    for (int i = 0; i < M; ++i) {
        if (i) cout << " ";
        cout << B[i];
    }
    cout << endl;
}

void Codeforces205Div1E() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    int N, M;
    cin >> N >> M;
    FPS<mint> C(M+1, 0);
    for (int i = 0; i < N; ++i) {
        int c;
        cin >> c;
        if (c > M) continue;
        C[c] += 1;
    }
    FPS<mint> F = inv(sqrt_base(C * mint(-4) + 1) + 1) * 2;
    for (int w = 1; w <= M; ++w) cout << F[w] << endl;
}

void TDPC_T() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;
    
    long long K, N;
    cin >> K >> N;
    
    --N;
    FPS<mint> P(K), Q(K + 1);
    Q[0] = 1;
    for (int i = 0; i < P.size(); ++i) P[i] = mint(1 - i);
    for (int i = 1; i < Q.size(); ++i) Q[i] = mint(-1);
    cout << BostanMori(P, Q, N) << endl;
}


int main() {
    //Yosupo_Inv_of_FPS();
    //Yosupo_Exp_of_FPS();
    //Yosupo_Log_of_FPS();
    Yosupo_Pow_of_FPS();
    //Yosupo_Sqrt_of_FPS();
    //HackerRankArrayRestoring();
    //Codeforces205Div1E();
    //TDPC_T();
}
