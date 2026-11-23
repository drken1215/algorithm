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

// int 128
using i128 = __int128_t;
i128 to_integer(const string &s) {
    i128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    if (s[0] == '-') res *= -1;
    return res;
}
istream& operator >> (istream &is, i128 &x) {
    string s;
    is >> s;
    x = to_integer(s);
    return is;
}
ostream& operator << (ostream &os, const i128 &x) {
    i128 ax = (x >= 0 ? x : -x);
    char buffer[128];
    char *d = end(buffer);
    do {
         --d;
        *d = "0123456789"[ax % 10];
        ax /= 10;
    } while (ax != 0);
    if (x < 0) {
        --d;
        *d = '-';
    }
    int len = end(buffer) - d;
    if (os.rdbuf()->sputn(d, len) != len) {
        os.setstate(ios_base::badbit);
    }
    return os;
}
i128 gcd(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    if (b == 0) return a;
    else return gcd(b, a % b);
}

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