//
// カタラン数 C(N) の計算
//　　・N × N グリッドにおいて、左下から右上まで、対角線の下のみを通る最短経路の本数
//
// verified:
//   KUPC 2019 D - Maximin Game
//     https://atcoder.jp/contests/kupc2019/tasks/kupc2019_d
//
//   SRM 401 DIV1 Easy
//     https://community.topcoder.com/stat?c=problem_statement&pm=8776&rd=12173
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

// Binomial coefficient (includes Catalan number)
template<class T> struct BiCoef {
    vector<T> fact_, inv_, finv_;
    constexpr BiCoef() {}
    constexpr BiCoef(int n) noexcept : fact_(n, 1), inv_(n, 1), finv_(n, 1) {
        init(n);
    }
    constexpr void init(int n) noexcept {
        fact_.assign(n, 1), inv_.assign(n, 1), finv_.assign(n, 1);
        int MOD = fact_[0].get_mod();
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
    
    // Catalan Number
    constexpr T catalan(int n) const noexcept {
        return com(n*2, n) - com(n*2, n-1);
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void KUPC_2019_D() {
    using mint = Fp<998244353>;
    int N;
    string S;
    cin >> N >> S;
    
    BiCoef<mint> bc(210000);
    mint res = 1;
    for (int i = 0; i < S.size();) {
        int j = i;
        while (j < S.size() && S[j] == S[i]) ++j;
        res *= bc.catalan(j - i);
        i = j;
    }
    cout << res << endl;
}

int main() {
    KUPC_2019_D();
}




