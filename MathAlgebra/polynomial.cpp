//
// 多項式アルゴリズム (by NTT, FPS)
// 　・係数は p を素数として Fp 体を想定
//
// verified:
//   Yosupo Library Checker - Convolution (Mod 1,000,000,007)
//     https://judge.yosupo.jp/problem/convolution_mod_1000000007
//
//   Yosupo Library Checker - Division of Polynomials
//     https://judge.yosupo.jp/problem/division_of_polynomials
//
//   Yosupo Library Checker - Polynomial Taylor Shift
//     https://judge.yosupo.jp/problem/polynomial_taylor_shift
//
//   AGC 005 F - Many Easy Problems (Polynomial Taylor Shift % mod 924844033)
//     https://atcoder.jp/contests/agc005/tasks/agc005_f
//
//   AOJ 2213 The Number of Solutions for a Polynomial (多項式の割り算 mod Any, 多項式 GCD)
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2213
//


#include <bits/stdc++.h>
using namespace std;


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

    // mul by convolution
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
        mint mod0 = MOD0, mod01 = mod0 * MOD1;
        vector<mint> res(N + M - 1);
        for (int i = 0; i < N + M - 1; ++i) {
            int y0 = c0[i].val;
            int y1 = (imod0 * (c1[i] - y0)).val;
            int y2 = (imod01 * (c2[i] - y0) - imod1 * y1).val;
            res[i] = mod01 * y2 + mod0 * y1 + y0;
        }
        return res;
    }
};

// Polynomial
template<typename mint> struct Poly : vector<mint> {
    using vector<mint>::vector;
 
    // constructor
    constexpr Poly(const vector<mint> &r) : vector<mint>(r) {}
 
    // core operator
    constexpr mint eval(const mint &v) {
        mint res = 0;
        for (int i = (int)this->size()-1; i >= 0; --i) {
            res *= v;
            res += (*this)[i];
        }
        return res;
    }
    constexpr Poly& normalize() {
        while (!this->empty() && this->back() == 0) this->pop_back();
        return *this;
    }
 
    // basic operator
    constexpr Poly operator - () const noexcept {
        Poly res = (*this);
        for (int i = 0; i < (int)res.size(); ++i) res[i] = -res[i];
        return res;
    }
    constexpr Poly operator + (const mint &v) const { return Poly(*this) += v; }
    constexpr Poly operator + (const Poly &r) const { return Poly(*this) += r; }
    constexpr Poly operator - (const mint &v) const { return Poly(*this) -= v; }
    constexpr Poly operator - (const Poly &r) const { return Poly(*this) -= r; }
    constexpr Poly operator * (const mint &v) const { return Poly(*this) *= v; }
    constexpr Poly operator * (const Poly &r) const { return Poly(*this) *= r; }
    constexpr Poly operator / (const mint &v) const { return Poly(*this) /= v; }
    constexpr Poly operator / (const Poly &r) const { return Poly(*this) /= r; }
    constexpr Poly operator % (const Poly &r) const { return Poly(*this) %= r; }
    constexpr Poly operator << (int x) const { return Poly(*this) <<= x; }
    constexpr Poly operator >> (int x) const { return Poly(*this) >>= x; }
    constexpr Poly& operator += (const mint &v) {
        if (this->empty()) this->resize(1);
        (*this)[0] += v;
        return *this;
    }
    constexpr Poly& operator += (const Poly &r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return this->normalize();
    }
    constexpr Poly& operator -= (const mint &v) {
        if (this->empty()) this->resize(1);
        (*this)[0] -= v;
        return *this;
    }
    constexpr Poly& operator -= (const Poly &r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return this->normalize();
    }
    constexpr Poly& operator *= (const mint &v) {
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= v;
        return *this;
    }
    constexpr Poly& operator *= (const Poly &r) {
        return *this = NTT::mul((*this), r);
    }
    constexpr Poly& operator <<= (int x) {
        Poly res(x, 0);
        res.insert(res.end(), begin(*this), end(*this));
        return *this = res;
    }
    constexpr Poly& operator >>= (int x) {
        Poly res;
        res.insert(res.end(), begin(*this) + x, end(*this));
        return *this = res;
    }
    
    // division, pow
    constexpr Poly& operator /= (const mint &v) {
        assert(v != 0);
        mint iv = modinv(v);
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= iv;
        return *this;
    }
    constexpr Poly& operator /= (const Poly &r) {
        assert(!r.empty());
        assert(r.back() != 0);
        this->normalize();
        if (this->size() < r.size()) {
            this->clear();
            return *this;
        }
        int need = (int)this->size() - (int)r.size() + 1;
        *this = (rev().pre(need) * r.rev().inner_inv(need)).pre(need).rev();
        return *this;
    }
    constexpr Poly& operator %= (const Poly &r) {
        assert(!r.empty());
        assert(r.back() != 0);
        this->normalize();
        Poly q = (*this) / r;
        return *this -= q * r;
    }
    
    // FPS functions
    constexpr Poly pre(int siz) const {
        return Poly(begin(*this), begin(*this) + min((int)this->size(), siz));
    }
    constexpr Poly rev() const {
        Poly res = *this;
        reverse(begin(res), end(res));
        return res;
    }
    
    // df/dx
    constexpr Poly diff() const {
        int n = (int)this->size();
        Poly res(n-1);
        for (int i = 1; i < n; ++i) res[i-1] = (*this)[i] * i;
        return res;
    }
    
    // \int f dx
    constexpr Poly integral() const {
        int n = (int)this->size();
        Poly res(n+1, 0);
        for (int i = 0; i < n; ++i) res[i+1] = (*this)[i] / (i+1);
        return res;
    }

    // inv(f), f[0] must not be 0
    constexpr Poly inner_inv(int deg) const {
        assert((*this)[0] != 0);
        if (deg < 0) deg = (int)this->size();
        Poly res({mint(1) / (*this)[0]});
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + res - res * res * pre(i << 1)).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    constexpr Poly inner_inv() const {
        return inner_inv((int)this->size());
    }
    
    // log(f) = \int f'/f dx, f[0] must be 1
    constexpr Poly inner_log(int deg) const {
        assert((*this)[0] == 1);
        Poly res = (diff() * inner_inv(deg)).integral();
        res.resize(deg);
        return res;
    }
    constexpr Poly inner_log() const {
        return inner_log((int)this->size());
    }
    
    // exp(f), f[0] must be 0
    constexpr Poly inner_exp(int deg) const {
        assert((*this)[0] == 0);
        Poly res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = res * (pre(i << 1) - res.inner_log(i << 1) + 1).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    constexpr Poly inner_exp() const {
        return inner_exp((int)this->size());
    }
    
    // pow(f) = exp(e * log f)
    constexpr Poly inner_pow(long long e, int deg) const {
        if (e == 0) {
            Poly res(deg, 0);
            res[0] = 1;
            return res;
        }
        long long i = 0;
        while (i < (int)this->size() && (*this)[i] == 0) ++i;
        if (i == (int)this->size() || i > (deg - 1) / e) return Poly(deg, 0);
        mint k = (*this)[i];
        Poly res = ((((*this) >> i) / k).inner_log(deg) * e).inner_exp(deg)
                    * mint(k).inner_pow(e) << (e * i);
        res.resize(deg);
        return res;
    }
    constexpr Poly inner_pow(long long e) const {
        return inner_pow(e, (int)this->size());
    }
};


/*/////////////////////////////*/
// Polynomial Algorithms
/*/////////////////////////////*/

// Binomial coefficient
template<class T> struct BiCoef {
    vector<T> fact_, inv_, finv_;
    constexpr BiCoef() {}
    constexpr BiCoef(int n) : fact_(n, 1), inv_(n, 1), finv_(n, 1) {
        init(n);
    }
    constexpr void init(int n) {
        fact_.assign(n, 1), inv_.assign(n, 1), finv_.assign(n, 1);
        int MOD = fact_[0].get_mod();
        for(int i = 2; i < n; i++){
            fact_[i] = fact_[i-1] * i;
            inv_[i] = -inv_[MOD%i] * (MOD/i);
            finv_[i] = finv_[i-1] * inv_[i];
        }
    }
    constexpr T com(int n, int k) const {
        if (n < k || n < 0 || k < 0) return 0;
        return fact_[n] * finv_[k] * finv_[n-k];
    }
    constexpr T fact(int n) const {
        if (n < 0) return 0;
        return fact_[n];
    }
    constexpr T inv(int n) const {
        if (n < 0) return 0;
        return inv_[n];
    }
    constexpr T finv(int n) const {
        if (n < 0) return 0;
        return finv_[n];
    }
};

// Polynomial Taylor Shift
// given: f(x), c
// find: coefficients of f(x + c)
template<class mint> Poly<mint> PolynomialTaylorShift(const Poly<mint> &f, long long c) {
    int N = (int)f.size() - 1;
    BiCoef<mint> bc(N + 1);
    
    // convolution
    Poly<mint> p(N + 1), q(N + 1);
    for (int i = 0; i <= N; ++i) {
        p[i] = f[i] * bc.fact(i);
        q[N - i] = mint(c).pow(i) * bc.finv(i);
    }
    Poly<mint> pq = p * q;
    
    // result
    Poly<mint> res(N + 1);
    for (int i = 0; i <= N; ++i) res[i] = pq[i + N] * bc.finv(i);
    return res;
}


/*/////////////////////////////*/
// for any mod
/*/////////////////////////////*/

// dynamic modint
struct DynamicModint {
    using mint = DynamicModint;
    
    // static menber
    static int MOD;
    
    // inner value
    long long val;
    
    // constructor
    DynamicModint() : val(0) { }
    DynamicModint(long long v) : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    long long get() const { return val; }
    static int get_mod() { return MOD; }
    static void set_mod(int mod) { MOD = mod; }
    
    // arithmetic operators
    mint operator + () const { return mint(*this); }
    mint operator - () const { return mint(0) - mint(*this); }
    mint operator + (const mint &r) const { return mint(*this) += r; }
    mint operator - (const mint &r) const { return mint(*this) -= r; }
    mint operator * (const mint &r) const { return mint(*this) *= r; }
    mint operator / (const mint &r) const { return mint(*this) /= r; }
    mint& operator += (const mint &r) {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    mint& operator -= (const mint &r) {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    mint& operator *= (const mint &r) {
        val = val * r.val % MOD;
        return *this;
    }
    mint& operator /= (const mint &r) {
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
    mint pow(long long n) const {
        mint res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    mint inv() const {
        mint res(1), div(*this);
        return res / div;
    }

    // other operators
    bool operator == (const mint &r) const {
        return this->val == r.val;
    }
    bool operator != (const mint &r) const {
        return this->val != r.val;
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
        is >> x.val;
        x.val %= x.get_mod();
        if (x.val < 0) x.val += x.get_mod();
        return is;
    }
    friend ostream& operator << (ostream &os, const mint &x) {
        return os << x.val;
    }
    friend mint pow(const mint &r, long long n) {
        return r.pow(n);
    }
    friend mint inv(const mint &r) {
        return r.inv();
    }
};

int DynamicModint::MOD;



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void Yosupo_Convolution_mod_1000000007() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;
    
    int N, M;
    cin >> N >> M;
    Poly<mint> a(N), b(M);
    for (int i = 0; i < N; ++i) cin >> a[i];
    for (int i = 0; i < M; ++i) cin >> b[i];
    Poly<mint> res = a * b;
    for (auto v : res) cout << v << " ";
    cout << endl;
}

void Yosupo_Division_of_Polynomials() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    int N, M;
    cin >> N >> M;
    Poly<mint> f(N), g(M);
    for (int i = 0; i < N; ++i) cin >> f[i];
    for (int i = 0; i < M; ++i) cin >> g[i];
    
    Poly<mint> q = f / g, r = f - g * q;
    cout << q.size() << " " << r.size() << endl;
    for (auto v : q) cout << v << " ";
    cout << endl;
    for (auto v : r) cout << v << " ";
    cout << endl;
}

void Yosupo_Polynomial_Taylor_Shift() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    long long N, c;
    cin >> N >> c;
    Poly<mint> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    
    Poly<mint> res = PolynomialTaylorShift(a, c);
    for (int i = 0; i < N; ++i) cout << res[i] << " ";
    cout << endl;
}

// AGC 005 F
void AGC_005_F() {
    const int MOD = 924844033;
    using mint = Fp<MOD>;
    
    int N;
    cin >> N;
    BiCoef<mint> bc(N + 1);
    vector<vector<int>> G(N);
    for (int i = 0; i < N-1; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        G[a].push_back(b);
        G[b].push_back(a);
    }
    
    Poly<mint> f(N+1, 0);
    vector<int> si(N, 1);
    auto rec = [&](auto self, int v, int p = -1) -> void {
        for (auto ch : G[v]) {
            if (ch == p) continue;
            self(self, ch, v);
            ++f[si[ch]];
            si[v] += si[ch];
        }
        ++f[N - si[v]];
    };
    rec(rec, 0);
    
    Poly<mint> g = PolynomialTaylorShift(f, 1);
    for (int k = 1; k <= N; ++k) cout << bc.com(N, k) * N - g[k] << endl;
}

// AOJ 2213
void AOJ_2213() {
    using mint = DynamicModint;
    
    auto pow = [&](const Poly<mint> &f, long long n, const Poly<mint> &m) -> Poly<mint> {
        Poly<mint> res = {1}, mul = f;
        while (n > 0) {
            if (n & 1) res = res * mul % m;
            mul = mul * mul % m;
            n >>= 1;
        }
        return res;
    };
    auto gcd = [&](auto self, const Poly<mint> &a, const Poly<mint> &b) -> Poly<mint> {
        if (b.empty()) return a;
        else return self(self, b, a % b);
    };
    
    long long N, P;
    while (cin >> N >> P) {
        if (N == 0) break;
        mint::set_mod(P);
        
        Poly<mint> f(N+1);
        for (int i = 0; i < N+1; ++i) cin >> f[i];
        f.normalize();
        if (f.empty()) {
            cout << P << endl;
            continue;
        }
        Poly<mint> x = {0, 1};
        Poly<mint> g = pow(x, P, f) - x;
        Poly<mint> res = gcd(gcd, f, g);
        cout << (int)res.size() - 1 << endl;
    }
}


int main() {
    //Yosupo_Convolution_mod_1000000007();
    //Yosupo_Division_of_Polynomials();
    //Yosupo_Polynomial_Taylor_Shift();
    //AGC_005_F();
    AOJ_2213();
}
