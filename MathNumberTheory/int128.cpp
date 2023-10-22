//
// 128 ビット整数 __int128 型のラッパー
//
// verified
//   Yosupo Library Checker - Sum of Floor of Linear
//     https://judge.yosupo.jp/problem/sum_of_floor_of_linear
//


#include <bits/stdc++.h>
using namespace std;


struct i128 {
    // inner value
    __int128 val;
    
    // constructor
    constexpr i128() {}
    constexpr i128(long long v) : val(v) {}
    constexpr i128(const string &s) {
        parse(s);
    }
    constexpr void parse(const string &s) {
        val = 0;
        for (char c : s) {
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


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void mini_test() {
    i128 a = 99;
    i128 b("993421434234324234234432421234324");
    
    cout << a + b << endl;
    cout << a * b << endl;
    cout << a % b << endl;
    cout << b * 2 << endl;
}

// Library Checker
// sum_{i=0}^{n-1} floor((a * i + b) / m)
template<class T> T floor_sum(T n, T a, T b, T m) {
    if (n == 0) return 0;
    T res = 0;
    if (a >= m) {
        res += n * (n - 1) * (a / m) / 2;
        a %= m;
    }
    if (b >= m) {
        res += n * (b / m);
        b %= m;
    }
    if (a == 0) return res;
    T ymax = (a * n + b) / m, xmax = ymax * m - b;
    if (ymax == 0) return res;
    res += (n - (xmax + a - 1) / a) * ymax;
    res += floor_sum(ymax, m, (a - xmax % a) % a, a);
    return res;
}

void YosupoSumOfFloorOfLinear() {
    int T;
    cin >> T;
    while (T--) {
        i128 N, M, A, B;
        cin >> N >> M >> A >> B;
        cout << floor_sum(N, A, B, M) << endl;
    }
}

int main() {
    //mini_test();
    YosupoSumOfFloorOfLinear();
}
