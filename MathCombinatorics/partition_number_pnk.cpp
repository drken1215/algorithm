//
// 分割数 P(n, k) を求める, O(nk) by 漸化式
//
// cf:
//  「写像12相」を総整理！ 〜 数え上げ問題の学びの宝庫 〜
//     https://qiita.com/drken/items/f2ea4b58b0d21621bd51
//
// verified:
//   yukicoder No.269 見栄っ張りの募金活動
//     https://yukicoder.me/problems/no/269
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

// partition number
template<class mint> struct PartionNumber {
    vector<vector<mint>> val;
    
    constexpr PartionNumber() {}
    constexpr PartionNumber(int MAX_N, int MAX_K) noexcept {
        init(MAX_N, MAX_K);
    }
    constexpr void init(int MAX_N, int MAX_K) noexcept {
        val.assign(MAX_N + 1, vector<mint>(MAX_K + 1, 0));
        val[0][0] = 1;
        for (int n = 0; n <= MAX_N; ++n) {
            for (int k = 1; k <= MAX_K; ++k) {
                val[n][k] += val[n][k - 1];
                if (n >= k) val[n][k] += val[n - k][k];
            }
        }
    }
    constexpr mint get(long long n, long long k) const noexcept {
        if (n < 0 || k <= 0 || n >= val.size() || k >= val[0].size()) return 0;
        return val[n][k];
    }
};

const int MOD = 1000000007;
using mint = Fp<MOD>;

int main() {
    long long N, S, K;
    cin >> N >> S >> K;
    
    PartionNumber<mint> pn(S, N);
    cout << pn.get(S - N * (N - 1) / 2 * K, N) << endl;
}

