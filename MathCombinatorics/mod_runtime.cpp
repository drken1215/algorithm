//
// 実行時に法が決まる Modular Arithmetics
//
// cf.
//   noshi91: modint 構造体を使ってみませんか？ (C++)
//     http://noshi91.hatenablog.com/entry/2019/03/31/174006
//
//   drken: 「1000000007 で割ったあまり」の求め方を総特集！ ~ 逆元から離散対数まで ~
//     https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a
//
// verified:
//   ARC 096 E - Everything on It
//     https://atcoder.jp/contests/arc096/tasks/arc096_c  
//


#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;


// modint
vector<int> MODS = { 1000000007 }; // 実行時に決まる
template<int IND = 0> struct Fp {
    long long val;
    
    int MOD = MODS[IND];
    constexpr Fp(long long v = 0) noexcept : val(v % MODS[IND]) {
        if (val < 0) val += MOD;
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
    friend constexpr ostream& operator << (ostream &os, const Fp<IND>& x) noexcept {
        return os << x.val;
    }
    friend constexpr istream& operator >> (istream &is, Fp<IND>& x) noexcept {
        return is >> x.val;
    }
    friend constexpr Fp<IND> modpow(const Fp<IND> &a, long long n) noexcept {
        if (n == 0) return 1;
        auto t = modpow(a, n / 2);
        t = t * t;
        if (n & 1) t = t * a;
        return t;
    }
};


// 二項係数ライブラリ
template<class T> struct BiCoef {
    vector<T> fact_, inv_, finv_;
    constexpr BiCoef(int n) noexcept : fact_(n, 1), inv_(n, 1), finv_(n, 1) {
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


// スターリング数 (n 個を k グループにわける、n >= k)
template<class T> struct Stirling {
    vector<vector<T> > S;
    constexpr Stirling(int MAX) noexcept : S(MAX, vector<T>(MAX, 0)) {
        S[0][0] = 1;
        for (int n = 1; n < MAX; ++n) {
            for (int k = 1; k <= n; ++k) {
                S[n][k] = S[n-1][k-1] + S[n-1][k] * k;
            }
        }
    }
    constexpr T get(int n, int k) {
        if (n < 0 || k < 0 || n < k) return 0;
        return S[n][k];
    }
};



const int MAX = 3100;
int main() {
    // 入力
    long long N;
    cin >> N >> MODS[0];
    using mint = Fp<>;
    
    // 前計算
    BiCoef<mint> bc(MAX); // 二項係数計算の前処理
    Stirling<mint> sl(MAX); // スターリング数の前処理

    // 2^n や 2^2^n の前計算、2^2^(n+1) = (2^2^n)^2
    vector<mint> two(MAX*MAX, 0), dtwo(MAX, 0);
    two[0] = 1, dtwo[0] = 2;
    for (int i = 1; i < MAX; ++i) dtwo[i] = dtwo[i-1] * dtwo[i-1];
    for (int i = 1; i < MAX*MAX; ++i) two[i] = two[i-1] * 2;
    
    // 求める
    mint res = 0;
    for (int n = 0; n <= N; ++n) {
        mint add = 0;
        for (int k = 0; k <= n; ++k) {
            mint jiyudo = two[(N-n)*k] * dtwo[N-n];
            mint core = sl.get(n, k) + sl.get(n, k+1) * (k+1);
            add += core * jiyudo;
        }
        mint choose = bc.com(N, n);
        add *= choose;
        if (n % 2 == 0) res += add;
        else res -= add;
    }
    cout << res << endl;
}
