//
// a^n mod. m
//
// cf.
//   drken: 「1000000007 で割ったあまり」の求め方を総特集！ 〜 逆元から離散対数まで 〜
//     https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a
//
// verified:
//   AOJ Course NTL_1_B Elementary Number Theory - Power
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=NTL_1_B&lang=jp
//


#include <iostream>
using namespace std;


long long modpow(long long a, long long n, long long mod) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}


int main() {
    const int MOD = 1000000007;
    long long m, n;
    cin >> m >> n;
    cout << modpow(m, n, MOD) << endl;
}
