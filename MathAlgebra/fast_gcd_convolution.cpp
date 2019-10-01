//
// Åº»ú GCD ¾ö¤ß¹þ¤ß
//
// verified:
//   yukicoder No.886 Direct
//     https://yukicoder.me/problems/no/886
// 


/*

  f(k) = sum_{GCD(i, j) = k} g(i)h(j)
  F(k) = sum_{k | i} f(i)

*/


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


template<int MOD> struct Fp {
    long long val;
    constexpr Fp(long long v = 0) noexcept : val(v % MOD) {
        if (val < 0) v += MOD;
    }
    constexpr int getmod() { return MOD; }
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
            a -= t * b; swap(a, b);
            u -= t * v; swap(u, v);
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
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD>& x) noexcept {
        return os << x.val;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD>& x) noexcept {
        return is >> x.val;
    }
    friend constexpr Fp<MOD> modpow(const Fp<MOD> &a, long long n) noexcept {
        if (n == 0) return 1;
        auto t = modpow(a, n / 2);
        t = t * t;
        if (n & 1) t = t * a;
        return t;
    }
};


// f(k) = sum_{GCD(i, j) = k} g(i)h(j)
// F(k) = sum_{k | i} f(i)
template<class T> struct FastGCDConvolution {
    int N;
    vector<bool> is_prime;

    FastGCDConvolution(int N) : N(N), is_prime(N, true) {
        is_prime[0] = is_prime[1] = false;
        for (int p = 2; p < N; ++p) {
            if (!is_prime[p]) continue;
            for (int i = p*2; i < N; i += p)
                is_prime[i] = false;
        }
    }

    void zeta(vector<T> &v) {
        for (int p = 2; p < N; ++p) {
            if (!is_prime[p]) continue;
            for (int i = (N-1)/p; i >= 1; --i)
                v[i] += v[i*p];
        }
    }

    void mebius(vector<T> &v) {
        for (int p = 2; p < N; ++p) {
            if (!is_prime[p]) continue;
            for (int i = 1; i*p < N; ++i)
                v[i] -= v[i*p];
        }
    }

    vector<T> mult(vector<T> &a, vector<T> &b) {
        vector<T> c(N);
        zeta(a); zeta(b); 
        for (int i = 1; i < N; ++i) c[i] = a[i] * b[i];
        mebius(c);
        return c;
    }
};


using mint = Fp<1000000007>;

int main() {
    int H, W; cin >> H >> W;
    mint res = mint(H) * (W-1) + mint(W) * (H-1);

    int N = max(H, W);
    vector<mint> h(N, 0), w(N, 0), f(N, 0);
    for (int i = 1; i < H; ++i) h[i] = H - i;
    for (int i = 1; i < W; ++i) w[i] = W - i;
    
    FastGCDConvolution<mint> fz(N);
    f = fz.mult(h, w);
    res += f[1] * 2;
    cout << res << endl;
}
