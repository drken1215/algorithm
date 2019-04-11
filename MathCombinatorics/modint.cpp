//
// Modular Arithmetics
//
// cf.
//   noshi91: modint 構造体を使ってみませんか？ (C++)
//     http://noshi91.hatenablog.com/entry/2019/03/31/174006
//
//   drken: 「1000000007 で割ったあまり」の求め方を総特集！ ~ 逆元から離散対数まで ~
//     https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a
//
// verified:
//   MUJIN 2018 F - チーム分け
//     https://atcoder.jp/contests/mujin-pc-2018/tasks/mujin_pc_2018_f  
// 


#include <iostream>
#include <vector>
using namespace std;


template<int MODULO> struct Fp {
    long long val;

    constexpr Fp(long long v = 0) noexcept : val(v % MODULO) {
        if (val < 0) v += MODULO;
    }
    constexpr Fp operator - () const noexcept {
        return val ? MODULO - val : 0;
    }
    constexpr Fp operator + (const Fp& r) const noexcept { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp& r) const noexcept { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp& r) const noexcept { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp& r) const noexcept { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp& r) noexcept {
        val += r.val;
        if (val >= MODULO) val -= MODULO;
        return *this;
    }
    constexpr Fp& operator -= (const Fp& r) noexcept {
        val -= r.val;
        if (val < 0) val += MODULO;
        return *this;
    }
    constexpr Fp& operator *= (const Fp& r) noexcept {
        val = val * r.val % MODULO;
        return *this;
    }
    constexpr Fp& operator /= (const Fp& r) noexcept {
        long long a = r.val, b = MODULO, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b; swap(a, b);
            u -= t * v; swap(u, v);
        }
        val = val * u % MODULO;
        if (val < 0) val += MODULO;
        return *this;
    }
    constexpr bool operator == (const Fp& r) const noexcept {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp& r) const noexcept {
        return this->val != r.val;
    }
};

template<int MODULO> constexpr ostream& operator <<
(ostream &os, const Fp<MODULO>& x) noexcept {
    return os << x.val;
}
template<int MODULO> constexpr istream& operator >>
(istream &is, Fp<MODULO>& x) noexcept {
    return is >> x.val;
}

template<int MODULO> constexpr Fp<MODULO> modpow
(const Fp<MODULO> &a, long long n) noexcept {
    if (n == 0) return 1;
    auto t = modpow(a, n / 2);
    t = t * t;
    if (n & 1) t = t * a;
    return t;
}


template<int MODULO> struct BiCoef {
    vector<Fp<MODULO> > fac, inv, finv;
    constexpr BiCoef(int n = 210000) noexcept : fac(n, 1), inv(n, 1), finv(n, 1) {
        for(int i = 2; i < n; i++){
            fac[i] = fac[i-1] * i;
            inv[i] = -inv[MODULO%i] * (MODULO/i);
            finv[i] = finv[i-1] * inv[i];
        }
    }
    constexpr Fp<MODULO> com(int n, int k) const noexcept {
        if (n < k || n < 0 || k < 0) return 0;
        return fac[n] * finv[k] * finv[n-k];
    }
};


const int MAX = 201010;
const int MOD = 998244353;
using mint = Fp<MOD>;

int main() {     
    BiCoef<MOD> bc(MAX);
    int N; cin >> N;
    vector<int> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];

    // nums[v] := v 人以上 OK な人数
    vector<long long> nums(N+2, 0);
    for (int i = 0; i < N; ++i) nums[a[i]]++;
    for (int i = N; i >= 0; --i) nums[i] += nums[i+1];

    // DP
    vector<vector<mint> > dp(N+2, vector<mint>(N+1, 0));
    dp[N+1][0] = 1;
    for (long long x = N; x >= 1; --x) {
        for (long long y = 0; y <= nums[x]; ++y) {
            for (long long k = 0; k <= N; ++k) {
                long long y2 = y - x * k;
                if (y2 < 0) break;
                if (y2 > nums[x+1]) continue;
                mint choose = bc.com(nums[x] - y2, x * k);
                mint fact = bc.fac[x*k] / modpow(bc.fac[x], k) * bc.finv[k];
                dp[x][y] += dp[x+1][y2] * choose * fact;
            }
        }
    }
    cout << dp[1][N] << endl;
}
