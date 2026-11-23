//
// 平方剰余 (Tonelli–Shanks 法)
//
// verified:
//   Yosupo Judge - Sqrt Mod
//     https://judge.yosupo.jp/problem/sqrt_mod
//


#include <bits/stdc++.h>
using namespace std;


// dynamic modint
struct DynamicModint {
    using mint = DynamicModint;
    
    // static menber
    static int MOD;
    
    // inner value
    unsigned int val;
    
    // constructor
    DynamicModint() : val(0) { }
    template<std::signed_integral T> DynamicModint(T v) {
        long long tmp = (long long)(v % (long long)(get_umod()));
        if (tmp < 0) tmp += get_umod();
        val = (unsigned int)(tmp);
    }
    template<std::unsigned_integral T> DynamicModint(T v) {
        val = (unsigned int)(v % get_umod());
    }
    long long get() const { return val; }
    static int get_mod() { return MOD; }
    static unsigned int get_umod() { return MOD; }
    static void set_mod(int mod) { MOD = mod; }
    
    // arithmetic operators
    mint operator + () const { return mint(*this); }
    mint operator - () const { return mint() - mint(*this); }
    mint operator + (const mint &r) const { return mint(*this) += r; }
    mint operator - (const mint &r) const { return mint(*this) -= r; }
    mint operator * (const mint &r) const { return mint(*this) *= r; }
    mint operator / (const mint &r) const { return mint(*this) /= r; }
    mint& operator += (const mint &r) {
        val += r.val;
        if (val >= get_umod()) val -= get_umod();
        return *this;
    }
    mint& operator -= (const mint &r) {
        val -= r.val;
        if (val >= get_umod()) val += get_umod();
        return *this;
    }
    mint& operator *= (const mint &r) {
        unsigned long long tmp = val;
        tmp *= r.val;
        val = (unsigned int)(tmp % get_umod());
        return *this;
    }
    mint& operator /= (const mint &r) {
        return *this = *this * r.inv(); 
    }
    mint pow(long long n) const {
        assert(n >= 0);
        mint res(1), mul(*this);
        while (n) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    mint inv() const {
        assert(val);
        assert(gcd(val, get_umod()) == 1);
        long long m = get_umod(), a = val, b = m, u = 1, v = 0;
        while (b > 0) {
            auto t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        return mint(u);
    }

    // other operators
    bool operator == (const mint &r) const {
        return this->val == r.val;
    }
    bool operator != (const mint &r) const {
        return this->val != r.val;
    }
    bool operator < (const mint &r) const {
        return this->val < r.val;
    }
    bool operator > (const mint &r) const {
        return this->val > r.val;
    }
    bool operator <= (const mint &r) const {
        return this->val <= r.val;
    }
    bool operator >= (const mint &r) const {
        return this->val >= r.val;
    }
    mint& operator ++ () {
        ++val;
        if (val == get_umod()) val = 0;
        return *this;
    }
    mint& operator -- () {
        if (val == 0) val = get_umod();
        --val;
        return *this;
    }
    mint operator ++ (int) {
        mint res = *this;
        ++*this;
        return res;
    }
    mint operator -- (int) {
        mint res = *this;
        --*this;
        return res;
    }
    friend istream& operator >> (istream &is, mint &x) {
        long long tmp = 1;
        is >> tmp;
        tmp = tmp % (long long)(get_umod());
        if (tmp < 0) tmp += get_umod();
        x.val = (unsigned int)(tmp);
        return is;
    }
    friend ostream& operator << (ostream &os, const mint &x) {
        return os << x.val;
    }
    friend mint pow(const mint &r, long long n) {
        return r.pow(n);
    }
    friend mint inv(const mint &r) {
        return r.inv();
    }
};
int DynamicModint::MOD;

// mod sqrt
template<class T_VAL, class T_MOD>
T_VAL mod_sqrt(T_VAL a, T_MOD p) {
    a = (a % p + p) % p;
    if (a <= 1) return a;
    using mint = DynamicModint;
    mint::set_mod(p);
    if (mint(a).pow((p - 1) >> 1) != 1) return T_VAL(-1);
    mint b = 1, one = 1;
    while (b.pow((p - 1) >> 1) == 1) b++;
    T_VAL m = p - 1, e = 0;
    while (m % 2 == 0) m >>= 1, e++;
    mint x = mint(a).pow((m - 1) >> 1);
    mint y = mint(a) * x * x;
    x *= a;
    mint z = mint(b).pow(m);
    while (y != 1) {
        T_VAL j = 0;
        mint t = y;
        while (t != one) {
            j++;
            t *= t;
        }
        z = z.pow(T_VAL(1) << (e - j - 1));
        x *= z, z *= z, y *= z;
        e = j;
    }
    T_VAL res = x.val;
    if (res * 2 > p) res = p - res;
    return res;
}


//------------------------------//
// Examples
//------------------------------//

// Library Checker - Sqrt MOD
void Yosupo_sqrt_mod() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int T;
    cin >> T;
    while (T--) {
        long long Y, P;
        cin >> Y >> P;
        long long res = mod_sqrt(Y, P);
        cout << res << endl;
    }
}


int main() {
    Yosupo_sqrt_mod();
}