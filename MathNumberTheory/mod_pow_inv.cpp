//
// mod_pow, mod_inv
//
// reference:
//   drken: 「998244353 で割ったあまり」の求め方を総特集！ ～ 逆元から離散対数まで ～
//     https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a
//
// verified:
//   simple test
//


#include <bits/stdc++.h>
using namespace std;


template<class T> T mod_pow(T a, T n, T m) {
    T res = 1;
    while (n > 0) {
        if (n % 2 == 1) res = res * a % m;
        a = a * a % m;
        n >>= 1;
    }
    return res;
};

template<class T> T mod_inv(T a, T m) {
    T b = m, u = 1, v = 0;
    while (b > 0) {
        T t = a / b;
        a -= t * b, swap(a, b);
        u -= t * v, swap(u, v);
    }
    u %= m;
    if (u < 0) u += m;
    return u;
};



//------------------------------//
// Examples
//------------------------------//

struct i128 {
    // inner value
    __int128 val;
    
    // constructor
    constexpr i128() : val(0) {}
    constexpr i128(long long v) : val(v) {}
    i128(const string &s) : val(0) {
        parse(s);
    }
    void parse(const string &s) {
        val = 0;
        for (auto c : s) {
            if (isdigit(c)) val = val * 10 + (c - '0');
        }
        if (s[0] == '-') val *= -1;
    }
    constexpr __int128 get() const {
        return val;
    }
    constexpr i128 abs() {
        if (val < 0) return -val;
        else return val;
    }
    
    // comparison operators
    constexpr bool operator == (const i128 &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const i128 &r) const {
        return this->val != r.val;
    }
    constexpr bool operator < (const i128 &r) const {
        return this->val < r.val;
    }
    constexpr bool operator > (const i128 &r) const {
        return this->val > r.val;
    }
    constexpr bool operator <= (const i128 &r) const {
        return this->val <= r.val;
    }
    constexpr bool operator >= (const i128 &r) const {
        return this->val >= r.val;
    }
    
    // arithmetic operators
    constexpr i128& operator += (const i128 &r) {
        val += r.val;
        return *this;
    }
    constexpr i128& operator -= (const i128 &r) {
        val -= r.val;
        return *this;
    }
    constexpr i128& operator *= (const i128 &r) {
        val *= r.val;
        return *this;
    }
    constexpr i128& operator /= (const i128 &r) {
        val /= r.val;
        return *this;
    }
    constexpr i128& operator %= (const i128 &r) {
        val %= r.val;
        return *this;
    }
    constexpr i128 operator + () const {
        return i128(*this);
    }
    constexpr i128 operator - () const {
        return i128(0) - i128(*this);
    }
    constexpr i128 operator + (const i128 &r) const {
        return i128(*this) += r;
    }
    constexpr i128 operator - (const i128 &r) const {
        return i128(*this) -= r;
    }
    constexpr i128 operator * (const i128 &r) const {
        return i128(*this) *= r;
    }
    constexpr i128 operator / (const i128 &r) const {
        return i128(*this) /= r;
    }
    constexpr i128 operator % (const i128 &r) const {
        return i128(*this) %= r;
    }
    
    // bit operators
    constexpr i128 operator >>= (long long r) {
        val >>= r;
        return *this;
    }
    constexpr i128 operator <<= (long long r) {
        val <<= r;
        return *this;
    }
    constexpr i128 operator &= (long long r) {
        val &= r;
        return *this;
    }
    constexpr i128 operator |= (long long r) {
        val |= r;
        return *this;
    }
    constexpr i128 operator << (long long r) const {
        return i128(*this) <<= r;
    }
    constexpr i128 operator >> (long long r) const {
        return i128(*this) >>= r;
    }
    constexpr i128 operator & (long long r) const {
        return i128(*this) &= r;
    }
    constexpr i128 operator | (long long r) const {
        return i128(*this) |= r;
    }
    
    // other operators
    constexpr i128& operator ++ () {
        ++val;
        return *this;
    }
    constexpr i128& operator -- () {
        --val;
        return *this;
    }
    constexpr i128 operator ++ (int) {
        i128 res = *this;
        ++*this;
        return res;
    }
    constexpr i128 operator -- (int) {
        i128 res = *this;
        --*this;
        return res;
    }
    friend istream& operator >> (istream &is, i128 &x) {
        string s;
        is >> s;
        x.parse(s);
        return is;
    }
    friend ostream& operator << (ostream &os, const i128 &x) {
        auto tmp = x.val < 0 ? -x.val : x.val;
        char buffer[128];
        char *d = end(buffer);
        do {
            --d;
            *d = "0123456789"[tmp % 10];
            tmp /= 10;
        } while (tmp != 0);
        if (x.val < 0) {
            --d;
            *d = '-';
        }
        int len = end(buffer) - d;
        if (os.rdbuf()->sputn(d, len) != len) {
            os.setstate(ios_base::badbit);
        }
        return os;
    }
};


int main() {
    const i128 P = 9999999999971LL;
    for (i128 a = 1; a <= 100; ++a) {
        cout << mod_pow(a, P - 2, P) << ", " << mod_inv(a, P) << endl;
    }
}

