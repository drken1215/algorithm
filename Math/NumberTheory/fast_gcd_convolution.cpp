//
// 添字 GCD 畳み込み
//
// verified:
//   AGC 038 C - LCMs
//     https://atcoder.jp/contests/agc038/tasks/agc038_c
//
//   yukicoder No.886 Direct
//     https://yukicoder.me/problems/no/886
// 


#include <bits/stdc++.h>
using namespace std;


// f(k) = sum_{GCD(i, j) = k} g(i)h(j)
// F(k) = sum_{k | i} f(i)
// F(k) = G(k)H(k)
namespace FastGCDConvolution {
    vector<bool> eratos(int N) {
        vector<bool> isprime(N, true);
        isprime[0] = isprime[1] = false;
        for (int p = 2; p < N; ++p) {
            if (!isprime[p]) continue;
            for (int i = p*2; i < N; i += p)
                isprime[i] = false;
        }
        return isprime;
    }

    template<class T> void zeta(vector<T> &v, const vector<bool> &isprime) {
        int N = (int)v.size();
        for (int p = 2; p < N; ++p) {
            if (!isprime[p]) continue;
            for (int i = (N-1)/p; i >= 1; --i)
                v[i] += v[i*p];
        }
    }

    template<class T> void mebius(vector<T> &v, const vector<bool> &isprime) {
        int N = (int)v.size();
        for (int p = 2; p < N; ++p) {
            if (!isprime[p]) continue;
            for (int i = 1; i*p < N; ++i)
                v[i] -= v[i*p];
        }
    }

    template<class T> vector<T> mul(const vector<T> &a, const vector<T> &b) {
        int N = max((int)a.size(), (int)b.size());
        const auto &isprime = eratos(N);
        vector<T> A(N, 0), B(N, 0);
        for (int i = 0; i < a.size(); ++i) A[i] = a[i];
        for (int i = 0; i < b.size(); ++i) B[i] = b[i];
        zeta(A, isprime), zeta(B, isprime);
        vector<T> C(N);
        for (int i = 1; i < N; ++i) C[i] = A[i] * B[i];
        mebius(C, isprime);
        return C;
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



//------------------------------//
// Examples
//------------------------------//

const int MOD = 998244353;
using mint = Fp<MOD>;

int main() {
    const int MAX = 1000001;
    int N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    vector<mint> F(MAX, 0);
    for (int i = 0; i < N; ++i) F[A[i]] += 1;
    for (int x = 1; x < MAX; ++x) F[x] *= x;
    auto FF = FastGCDConvolution::mul(F, F);
    mint res = 0;
    for (int x = 1; x < MAX; ++x) res += FF[x] / x;
    for (int i = 0; i < N; ++i) res -= A[i];
    cout << res / 2 << endl;
}
