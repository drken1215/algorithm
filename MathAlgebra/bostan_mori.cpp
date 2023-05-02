//
// Bostan-Mori 法
//   find [x^N] P(x)/Q(x), where deg(Q(x)) = K, deg(P(x)) < K
//   time complexity: O(K log K log N)
//
// verified:
//   TDPC T - フィボナッチ
//     https://atcoder.jp/contests/tdpc/tasks/tdpc_fibonacci
//


#include <bits/stdc++.h>
using namespace std;


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
    friend constexpr istream& operator >> (istream& is, Fp<MOD>& x) noexcept {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream& os, const Fp<MOD>& x) noexcept {
        return os << x.val;
    }
    friend constexpr Fp<MOD> modpow(const Fp<MOD>& r, long long n) noexcept {
        if (n == 0) return 1;
        if (n < 0) return modpow(modinv(r), -n);
        auto t = modpow(r, n / 2);
        t = t * t;
        if (n & 1) t = t * r;
        return t;
    }
    friend constexpr Fp<MOD> modinv(const Fp<MOD>& r) noexcept {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        return Fp<MOD>(u);
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
    template<class mint> void trans(vector<mint>& v, bool inv = false) {
        if (v.empty()) return;
        int N = (int)v.size();
        int MOD = v[0].getmod();
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
        if (min(N, M) < 30) return naive_mul(A, B);
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
    vector<long long> mul_ll
    (const vector<long long>& A, const vector<long long>& B) {
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
        FPS res = exp(log((f >> i) / k, deg) * e, deg) * modpow(k, e) << (e * i);
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

// Bostan-Mori
// find [x^N] P(x)/Q(x), O(K log K log N)
// deg(Q(x)) = K, deg(P(x)) < K
template <typename mint> mint BostanMori(const FPS<mint> &P, const FPS<mint> &Q, long long N) {
    assert(!P.empty() && !Q.empty());
    if (N == 0) return P[0] / Q[0];
    
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


////////////////////////////////////////
// solver
////////////////////////////////////////

void TDPC_T() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;
    
    // 入力
    long long K, N;
    cin >> K >> N;
    
    // Bostan-Mori
    FPS<mint> P(K), Q(K + 1);
    Q[0] = 1;
    for (int i = 0; i < P.size(); ++i) P[i] = mint(1 - i);
    for (int i = 1; i < Q.size(); ++i) Q[i] = mint(-1);
    cout << BostanMori(P, Q, N - 1) << endl;
}

int main() {
    TDPC_T();
}
