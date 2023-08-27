//
// 分割数 P(n) を O(N√N) で求める
//
// cf.
//   drken: 「写像12相」を総整理！ 〜 数え上げ問題の学びの宝庫 〜
//     https://qiita.com/drken/items/f2ea4b58b0d21621bd51
//
// verified:
//   Yosupo Library Checker - Partition Function
//     https://judge.yosupo.jp/problem/partition_function
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

// partition number
template<class T> struct PartitionFunction {
    vector<T> P;
    constexpr PartitionFunction(int MAX) noexcept : P(MAX+1, 0) {
        P[0] = 1;
        for (int i = 1; i <= MAX; ++i) {
            for (int j = 1, sign = 1; i-(j*j*3-j)/2 >= 0; ++j, sign *= -1) {
                P[i] += P[i-(j*j*3-j)/2] * sign;
                if (i-(j*j*3+j)/2 >= 0) P[i] += P[i-(j*j*3+j)/2] * sign;
            }
        }
    }
    constexpr T get(int n) {
        if (n < 0) return 0;
        return P[n];
    }
};


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void Yosupo_Partition_Function() {
    using mint = Fp<998244353>;
    int N;
    cin >> N;
    PartitionFunction<mint> pfun(N);
    for (int i = 0; i <= N; ++i) cout << pfun.get(i) << " ";
    cout << endl;
}

int main() {
    Yosupo_Partition_Function();
}


