//
// Formal Power Series (on runtime mod)
//
// verified:
//   Yosupo Judge
//     https://judge.yosupo.jp/problem/inv_of_formal_power_series
//
//   AOJ 2213 The Number of Solutions for a Polynomial
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2213
//


#include <bits/stdc++.h>
using namespace std;

#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, vector<vector<T> > P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, set<T> P)
{ for(auto it : P) { s << "<" << it << "> "; } return s << endl; }
template<class T1, class T2> ostream& operator << (ostream &s, map<T1,T2> P)
{ for(auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s << endl; }


// modint (replace MODS[0] on runtime)
vector<int> MODS = { 1000000007, 754974721, 167772161, 469762049 };
template<int IND = 0> struct Fp {
    long long val;
    constexpr Fp(long long v = 0) noexcept : val(v % MODS[IND]) {
        if (val < 0) val += MODS[IND];
    }
    constexpr int getmod() const { return MODS[IND]; }
    constexpr Fp operator - () const noexcept {
        return val ? MODS[IND] - val : 0;
    }
    constexpr Fp operator + (const Fp& r) const noexcept { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp& r) const noexcept { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp& r) const noexcept { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp& r) const noexcept { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp& r) noexcept {
        val += r.val;
        if (val >= MODS[IND]) val -= MODS[IND];
        return *this;
    }
    constexpr Fp& operator -= (const Fp& r) noexcept {
        val -= r.val;
        if (val < 0) val += MODS[IND];
        return *this;
    }
    constexpr Fp& operator *= (const Fp& r) noexcept {
        val = val * r.val % MODS[IND];
        return *this;
    }
    constexpr Fp& operator /= (const Fp& r) noexcept {
        long long a = r.val, b = MODS[IND], u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b; swap(a, b);
            u -= t * v; swap(u, v);
        }
        val = val * u % MODS[IND];
        if (val < 0) val += MODS[IND];
        return *this;
    }
    constexpr bool operator == (const Fp& r) const noexcept {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp& r) const noexcept {
        return this->val != r.val;
    }
    friend constexpr istream& operator >> (istream& is, Fp<IND>& x) noexcept {
        is >> x.val;
        x.val %= MODS[IND];
        if (x.val < 0) x.val += MODS[IND];
        return is;
    }
    friend constexpr ostream& operator << (ostream& os, const Fp<IND>& x) noexcept {
        return os << x.val;
    }
    friend constexpr Fp<IND> modpow(const Fp<IND>& r, long long n) noexcept {
        if (n == 0) return 1;
        auto t = modpow(r, n / 2);
        t = t * t;
        if (n & 1) t = t * r;
        return t;
    }
    friend constexpr Fp<IND> modinv(const Fp<IND>& r) noexcept {
        long long a = r.val, b = MODS[IND], u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        return Fp<IND>(u);
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
    set<int> setmod;
    template<class mint> void trans(vector<mint>& v, bool inv = false) {
        if (v.empty()) return;
        int N = (int)v.size();
        int MOD = v[0].getmod();
        int PR = calc_primitive_root(MOD);
        static vector<long long> vbw(30), vibw(30);
        if (!setmod.count(MOD)) {
            setmod.insert(MOD);
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
    static constexpr int IND1 = 1; //754974721;
    static constexpr int IND2 = 2; //167772161;
    static constexpr int IND3 = 3; //469762049;
    using mint1 = Fp<IND1>;
    using mint2 = Fp<IND2>;
    using mint3 = Fp<IND3>;
    static const mint2 imod1 = 95869806; // modinv(MOD1, MOD2);
    static const mint3 imod2 = 104391568; // modinv(MOD2, MOD3);
    static const mint3 imod12 = 187290749; // imod2 / MOD1;

    // small case (T = mint, long long)
    template<class T> vector<T> naive_mul 
    (const vector<T>& A, const vector<T>& B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        vector<T> res(N + M - 1);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < M; ++j)
                res[i + j] += A[i] * B[j];
        return res;
    }

    // mint
    template<class mint> vector<mint> mul
    (const vector<mint>& A, const vector<mint>& B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        if (min(N, M) < 80) return naive_mul(A, B);
        int MOD = A[0].getmod();
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
        vector<mint1> a1(size_fft, 0), b1(size_fft, 0), c1(size_fft, 0);
        vector<mint2> a2(size_fft, 0), b2(size_fft, 0), c2(size_fft, 0);
        vector<mint3> a3(size_fft, 0), b3(size_fft, 0), c3(size_fft, 0);
        for (int i = 0; i < N; ++i)
            a1[i] = A[i].val, a2[i] = A[i].val, a3[i] = A[i].val;
        for (int i = 0; i < M; ++i)
            b1[i] = B[i].val, b2[i] = B[i].val, b3[i] = B[i].val;
        trans(a1), trans(a2), trans(a3), trans(b1), trans(b2), trans(b3);
        for (int i = 0; i < size_fft; ++i) {
            c1[i] = a1[i] * b1[i];
            c2[i] = a2[i] * b2[i];
            c3[i] = a3[i] * b3[i];
        }
        trans(c1, true), trans(c2, true), trans(c3, true);
        mint mod1 = MODS[IND1], mod12 = mod1 * MODS[IND2];
        vector<mint> res(N + M - 1);
        for (int i = 0; i < N + M - 1; ++i) {
            int y1 = c1[i].val;
            int y2 = (imod1 * (c2[i] - y1)).val;
            int y3 = (imod12 * (c3[i] - y1) - imod2 * y2).val;
            res[i] = mod12 * y3 + mod1 * y2 + y1;
        }
        return res;
    }
};

// Binomial Coefficient
template<class T> struct BiCoef {
    vector<T> fact_, inv_, finv_;
    constexpr BiCoef() {}
    constexpr BiCoef(int n) noexcept : fact_(n, 1), inv_(n, 1), finv_(n, 1) {
        init(n);
    }
    constexpr void init(int n) noexcept {
        fact_.assign(n, 1), inv_.assign(n, 1), finv_.assign(n, 1);
        int MOD = fact_[0].getmod();
        for(int i = 2; i < n; i++){
            fact_[i] = fact_[i-1] * i;
            inv_[i] = -inv_[MOD%i] * (MOD/i);
            finv_[i] = finv_[i-1] * inv_[i];
        }
    }
    constexpr T com(int n, int k) const noexcept {
        if (n < k || n < 0 || k < 0) return 0;
        return fact_[n] * finv_[k] * finv_[n-k];
    }
    constexpr T fact(int n) const noexcept {
        if (n < 0) return 0;
        return fact_[n];
    }
    constexpr T inv(int n) const noexcept {
        if (n < 0) return 0;
        return inv_[n];
    }
    constexpr T finv(int n) const noexcept {
        if (n < 0) return 0;
        return finv_[n];
    }
};

// Formal Power Series
template <typename mint> struct FPS : vector<mint> {
    using vector<mint>::vector;
 
    // constructor
    FPS(const vector<mint>& r) : vector<mint>(r) {}
 
    // core operator
    inline FPS pre(int siz) const {
        return FPS(begin(*this), begin(*this) + min((int)this->size(), siz));
    }
    inline FPS rev() const {
        FPS res = *this;
        reverse(begin(res), end(res));
        return res;
    }
    inline FPS& normalize() {
        while (!this->empty() && this->back() == 0) this->pop_back();
        return *this;
    }
 
    // basic operator
    inline FPS operator - () const noexcept {
        FPS res = (*this);
        for (int i = 0; i < (int)res.size(); ++i) res[i] = -res[i];
        return res;
    }
    inline FPS operator + (const mint& v) const { return FPS(*this) += v; }
    inline FPS operator + (const FPS& r) const { return FPS(*this) += r; }
    inline FPS operator - (const mint& v) const { return FPS(*this) -= v; }
    inline FPS operator - (const FPS& r) const { return FPS(*this) -= r; }
    inline FPS operator * (const mint& v) const { return FPS(*this) *= v; }
    inline FPS operator * (const FPS& r) const { return FPS(*this) *= r; }
    inline FPS operator / (const mint& v) const { return FPS(*this) /= v; }
    inline FPS operator << (int x) const { return FPS(*this) <<= x; }
    inline FPS operator >> (int x) const { return FPS(*this) >>= x; }
    inline FPS& operator += (const mint& v) {
        if (this->empty()) this->resize(1);
        (*this)[0] += v;
        return *this;
    }
    inline FPS& operator += (const FPS& r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return this->normalize();
    }
    inline FPS& operator -= (const mint& v) {
        if (this->empty()) this->resize(1);
        (*this)[0] -= v;
        return *this;
    }
    inline FPS& operator -= (const FPS& r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return this->normalize();
    }
    inline FPS& operator *= (const mint& v) {
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= v;
        return *this;
    }
    inline FPS& operator *= (const FPS& r) {
        return *this = NTT::mul((*this), r);
    }
    inline FPS& operator /= (const mint& v) {
        assert(v != 0);
        mint iv = modinv(v);
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= iv;
        return *this;
    }
    inline FPS& operator <<= (int x) {
        FPS res(x, 0);
        res.insert(res.end(), begin(*this), end(*this));
        return *this = res;
    }
    inline FPS& operator >>= (int x) {
        FPS res;
        res.insert(res.end(), begin(*this) + x, end(*this));
        return *this = res;
    }
    inline mint eval(const mint& v){
        mint res = 0;
        for (int i = (int)this->size()-1; i >= 0; --i) {
            res *= v;
            res += (*this)[i];
        }
        return res;
    }
    inline friend FPS gcd(const FPS& f, const FPS& g) {
        if (g.empty()) return f;
        return gcd(g, f % g);
    }

    // advanced operation
    // df/dx
    inline friend FPS diff(const FPS& f) {
        int n = (int)f.size();
        FPS res(n-1);
        for (int i = 1; i < n; ++i) res[i-1] = f[i] * i;
        return res;
    }

    // \int f dx
    inline friend FPS integral(const FPS& f) {
        int n = (int)f.size();
        FPS res(n+1, 0);
        for (int i = 0; i < n; ++i) res[i+1] = f[i] / (i+1);
        return res;
    }

    // inv(f), f[0] must not be 0
    inline friend FPS inv(const FPS& f, int deg) {
        assert(f[0] != 0);
        if (deg < 0) deg = (int)f.size();
        FPS res({mint(1) / f[0]});
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + res - res * res * f.pre(i << 1)).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS inv(const FPS& f) {
        return inv(f, f.size());
    }

    // division, r must be normalized (r.back() must not be 0)
    inline FPS& operator /= (const FPS& r) {
        assert(!r.empty());
        assert(r.back() != 0);
        this->normalize();
        if (this->size() < r.size()) {
            this->clear();
            return *this;
        }
        int need = (int)this->size() - (int)r.size() + 1;
        *this = ((*this).rev().pre(need) * inv(r.rev(), need)).pre(need).rev();
        return *this;
    }
    inline FPS& operator %= (const FPS &r) {
        assert(!r.empty());
        assert(r.back() != 0);
        this->normalize();
        FPS q = (*this) / r;
        return *this -= q * r;
    }
    inline FPS operator / (const FPS& r) const { return FPS(*this) /= r; }
    inline FPS operator % (const FPS& r) const { return FPS(*this) %= r; }

    // log(f) = \int f'/f dx, f[0] must be 1
    inline friend FPS log(const FPS& f, int deg) {
        assert(f[0] == 1);
        FPS res = integral(diff(f) * inv(f, deg));
        res.resize(deg);
        return res;
    }
    inline friend FPS log(const FPS& f) {
        return log(f, f.size());
    }

    // exp(f), f[0] must be 0
    inline friend FPS exp(const FPS& f, int deg) {
        assert(f[0] == 0);
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = res * (f.pre(i<<1) - log(res, i<<1) + 1).pre(i<<1);
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS exp(const FPS& f) {
        return exp(f, f.size());
    }

    // pow(f) = exp(e * log f)
    inline friend FPS pow(const FPS& f, long long e, int deg) {
        long long i = 0;
        while (i < (int)f.size() && f[i] == 0) ++i;
        if (i == (int)f.size()) return FPS(deg, 0);
        if (i * e >= deg) return FPS(deg, 0);
        mint k = f[i];
        FPS res = exp(log((f >> i) / k) * e) * modpow(k, e) << (e * i);
        res.resize(deg);
        return res;
    }
    inline friend FPS pow(const FPS& f, long long e) {
        return pow(f, e, f.size());
    }

    // sqrt(f), f[0] must be 1
    inline friend FPS sqrt_base(const FPS& f, int deg) {
        assert(f[0] == 1);
        mint inv2 = mint(1) / 2;
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + f.pre(i << 1) * inv(res, i << 1)).pre(i << 1);
            for (mint& x : res) x *= inv2;
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS sqrt_base(const FPS& f) {
        return sqrt_base(f, f.size());
    }
};

using mint = Fp<>;



////////////////////////////////////////
// solver
////////////////////////////////////////

FPS<mint> modpow(const FPS<mint> &f, long long n, const FPS<mint> &m) {
    if (n == 0) return FPS<mint>(1, 1);
    auto t = modpow(f, n / 2, m);
    t = (t * t) % m;
    if (n & 1) t = (t * f) % m;
    auto q = t / m;
    auto r = t % m;
    return t;
}

void AOJ2213() {
    long long N, P;
    while (cin >> N >> P) {
        if (N == 0) break;
        MODS[0] = P;
		
        FPS<mint> f(N+1);
        for (int i = 0; i < N+1; ++i) cin >> f[i];
        f.normalize();
        if (f.empty()) { 
            cout << P << endl;
            continue; 
        }
        
        FPS<mint> x(2, 0);
        x[1] = 1; // ¿�༰ x
        FPS<mint> g = modpow(x, P, f) - x;
        FPS<mint> res = gcd(f, g);
        cout << (int)res.size() - 1 << endl;
    }
}

int main() {
    AOJ2213();
}
