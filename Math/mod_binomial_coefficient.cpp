//
// nCk mod p (オーソドックス ver.)
//
// cf.
//   drken: 「998244353 で割ったあまり」の求め方を総特集！ 〜 逆元から離散対数まで 〜
//     https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a
//
// verified:
//   ABC 127 E - Cell Distance
//     https://atcoder.jp/contests/abc127/tasks/abc127_e
//
//   AGC 025 B - RGB Coloring
//     https://atcoder.jp/contests/agc025/tasks/agc025_b
//


#include <bits/stdc++.h>
using namespace std;


// modint
template<int MOD = 998244353, bool PRIME = true> struct Fp {
    // inner value
    unsigned int val;
    
    // constructor
    constexpr Fp() : val(0) { }
    template<std::signed_integral T> constexpr Fp(T v) {
        long long tmp = (long long)(v % (long long)(get_umod()));
        if (tmp < 0) tmp += get_umod();
        val = (unsigned int)(tmp);
    }
    template<std::unsigned_integral T> constexpr Fp(T v) {
        val = (unsigned int)(v % get_umod());
    }
    constexpr long long get() const { return val; }
    constexpr static int get_mod() { return MOD; }
    constexpr static unsigned int get_umod() { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator + () const { return Fp(*this); }
    constexpr Fp operator - () const { return Fp() - Fp(*this); }
    constexpr Fp operator + (const Fp &r) const { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp &r) {
        val += r.val;
        if (val >= get_umod()) val -= get_umod();
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) {
        val -= r.val;
        if (val >= get_umod()) val += get_umod();
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) {
        unsigned long long tmp = val;
        tmp *= r.val;
        val = (unsigned int)(tmp % get_umod());
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) {
        return *this = *this * r.inv(); 
    }
    constexpr Fp pow(long long n) const {
        assert(n >= 0);
        Fp res(1), mul(*this);
        while (n) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const {
        assert(val);
        if (PRIME) {
            return pow(get_umod() - 2);
        } else {
            assert(gcd(val, get_umod()) == 1);
            long long m = get_umod(), a = val, b = m, u = 1, v = 0;
            while (b > 0) {
                auto t = a / b;
                a -= t * b, swap(a, b);
                u -= t * v, swap(u, v);
            }
            return Fp(u);
        }
    }

    // other operators
    constexpr bool operator == (const Fp &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const {
        return this->val != r.val;
    }
    constexpr bool operator < (const Fp &r) const {
        return this->val < r.val;
    }
    constexpr bool operator > (const Fp &r) const {
        return this->val > r.val;
    }
    constexpr bool operator <= (const Fp &r) const {
        return this->val <= r.val;
    }
    constexpr bool operator >= (const Fp &r) const {
        return this->val >= r.val;
    }
    constexpr Fp& operator ++ () {
        ++val;
        if (val == get_umod()) val = 0;
        return *this;
    }
    constexpr Fp& operator -- () {
        if (val == 0) val = get_umod();
        --val;
        return *this;
    }
    constexpr Fp operator ++ (int) {
        Fp res = *this;
        ++*this;
        return res;
    }
    constexpr Fp operator -- (int) {
        Fp res = *this;
        --*this;
        return res;
    }
    friend constexpr istream& operator >> (istream &is, Fp &x) {
        long long tmp = 1;
        is >> tmp;
        tmp = tmp % (long long)(get_umod());
        if (tmp < 0) tmp += get_umod();
        x.val = (unsigned int)(tmp);
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp &x) {
        return os << x.val;
    }
    friend constexpr Fp pow(const Fp &r, long long n) {
        return r.pow(n);
    }
    friend constexpr Fp inv(const Fp &r) {
        return r.inv();
    }
};

// Binomial coefficient
template<class mint> struct BiCoef {
    vector<mint> fact_, inv_, finv_;
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
    constexpr mint com(int n, int k) const {
        if (n < k || n < 0 || k < 0) return 0;
        return fact_[n] * finv_[k] * finv_[n-k];
    }
    constexpr mint fact(int n) const {
        if (n < 0) return 0;
        return fact_[n];
    }
    constexpr mint inv(int n) const {
        if (n < 0) return 0;
        return inv_[n];
    }
    constexpr mint finv(int n) const {
        if (n < 0) return 0;
        return finv_[n];
    }
};


//------------------------------//
// Examples
//------------------------------//

// ABC 127 E - Cell Distance
void ABC_127_E() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;

    long long N, M, K;
    cin >> N >> M >> K;

    BiCoef<mint> bc(N * M);
    mint sum = 0;
    for (int i = 0; i <= N-1; ++i) {
        for (int j = 0; j <= M-1; ++j) {
            mint tmp = mint(N - i) * mint(M - j) * mint(i + j);
            if (i != 0 && j != 0) tmp *= 2;
            sum += tmp;
        }
    }
    cout << sum * bc.com(N * M - 2, K - 2) << endl;
}

// AGC 025 B - RGB Coloring
void AGC_025_B() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;

    long long N, A, B, K;
    cin >> N >> A >> B >> K;

    BiCoef<mint> bc(N + 1);
    mint res = 0;
    for (long long a = 0; a <= N; ++a) { // A の個数を全探索
        long long rem = K - a * A;
        if (rem % B != 0) continue;
        long long b = rem / B;
        if (b > N) continue;
        mint tmp = bc.com(N, a) * bc.com(N, b);
        res += tmp;
    }   
    cout << res << endl;
}


int main() {
    //ABC_127_E();
    AGC_025_B();
}