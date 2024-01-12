//
// extended GCD (Euclid's algorithm)
//
// cf.
//   拡張ユークリッドの互除法 〜 一次不定方程式 ax + by = c の解き方 〜
//     https://qiita.com/drken/items/b97ff231e43bce50199a
//
// verified
//   AOJ Course NTL_1_E Elementary Number Theory - Extended Euclid Algorithm
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=NTL_1_E&lang=jp
//


#include <iostream>
using namespace std;


// 返り値: a と b の最大公約数
// ax + by = gcd(a, b) を満たす (x, y) が格納される
long long extGCD(long long a, long long b, long long &x, long long &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    long long d = extGCD(b, a%b, y, x);
    y -= a/b * x;
    return d;
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    long long a, b; cin >> a >> b;
    long long x, y;
    extGCD(a, b, x, y);
    cout << x << " " << y << endl;
}
