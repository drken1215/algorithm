//
// extended GCD (Euclid's algorithm)
//
// cf.
//   拡張ユークリッドの互除法 〜 一次不定方程式 ax + by = c の解き方 〜
//     https://qiita.com/drken/items/b97ff231e43bce50199a
//
// verified
//   AtCoder ABC 340 F - S = 1
//     https://atcoder.jp/contests/abc340/tasks/abc340_f
//
//   AOJ Course NTL_1_E Elementary Number Theory - Extended Euclid Algorithm
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=NTL_1_E&lang=jp
//


#include <iostream>
using namespace std;


// 返り値: a と b の最大公約数
// ax + by = gcd(a, b) を満たす (x, y) が格納される
template<class T> T ext_gcd(T a, T b, T &x, T &y) {
    T xsign = 1, ysign = 1;
    if (a < 0) {
        a = -a, xsign *= -1;
    }
    if (b < 0) {
        b = -b, ysign *= -1;
    }
    
    if (b == 0) {
        x = xsign, y = 0;
        return a;
    }
    T d = ext_gcd(b, a % b, y, x);
    y -= a / b * x;
    x *= xsign, y *= ysign;
    return d;
}



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

void ABC_340_F() {
    i128 X, Y;
    cin >> X >> Y;
    
    i128 a, b;
    i128 g = ext_gcd(X, Y, a, b);
    if (g == 1) {
        cout << b * 2 << " " << a * (-2) << endl;
    } else if (g == 2) {
        cout << b << " " << -a << endl;
    } else {
        cout << -1 << endl;
    }
}


int main() {
    ABC_340_F();
}


