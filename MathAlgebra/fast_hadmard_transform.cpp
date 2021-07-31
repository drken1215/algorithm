//
// 添字 XOR 畳み込み (高速アダマール変換)
//
// verified:
//   ABC 212 H - Nim Counting 
//     https://atcoder.jp/contests/abc212/tasks/abc212_h
//


#include <bits/stdc++.h>
using namespace std;


// Fast Hadamard Transform
// N must be 2^K for some K
namespace FastHadamardTransform {
    template<class T> void trans(vector<T> &v, bool inv = false) {
        int N = v.size();
        for (int i = 1; i < N; i <<= 1) {
            for (int j = 0; j < N; ++j) {
                if ((j & i) == 0) {
                    auto x = v[j], y = v[j | i];
                    v[j] = x + y;
                    v[j | i] = x - y;
                    if (inv) v[j] /= 2, v[j | i] /= 2;
                }
            }
        }
    }

    template<class T> vector<T> mul(const vector<T> &a, const vector<T> &b) {
        int N = a.size();
        auto A = a, B = b;
        trans(A), trans(B);
        vector<T> C(N);
        for (int i = 0; i < N; ++i) C[i] = A[i] * B[i];
        trans(C, true);
        return C;
    }
};
using namespace FastHadamardTransform;


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

const int MOD = 998244353;
using mint = Fp<MOD>;


// 級数和を求める
vector<mint> operator + (const vector<mint> &f, const vector<mint> &g) {
    vector<mint> res(f.size());
    for (int i = 0; i < f.size(); ++i) {
        res[i] = f[i] + g[i];
    }
    return res;
} 
vector<mint> operator * (const vector<mint> &f, const vector<mint> &g) {
    vector<mint> res(f.size());
    for (int i = 0; i < f.size(); ++i) {
        res[i] = f[i] * g[i];
    }
    return res;
}

// 級数和
// unit: 冪乗演算の単位元
template<class T> T calc_series(T v, long long N, T unit) {
    if (N == 1) return unit;

    if (N % 2 == 1)
        return v * calc_series(v, N - 1, unit) + unit;
    else
        return (v + unit) * calc_series(v * v, N / 2, unit);
}    

// ベクトルサイズ
const int V = 65536;

int main() {
    // 入力
    int N, K;
    cin >> N >> K;
    vector<mint> f(V);
    for (int i = 0; i < K; ++i) {
        int a;
        cin >> a;
        f[a] = 1;
    }

    // 全体の答え
    mint all = mint(K) * calc_series(mint(K), N, mint(1));

    // 級数和を求める
    vector<mint> unit(V, 0);
    unit[0] = 1;
    trans(f), trans(unit);  // 高速アダマール変換
    auto res = f * calc_series(f, N, unit);
    trans(res, true);

    cout << all - res[0] << endl;
}
