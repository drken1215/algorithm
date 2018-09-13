//
// べき乗
//
// 参考:
//   https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a#4-%E7%B4%AF%E4%B9%97-an
//
// verified
//   AOJ Course NTL_1_B Elementary Number Theory - Power
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=NTL_1_B&lang=jp
//

/*
    a^n mod. mod を計算する、O(log n)
 */




#include <iostream>
using namespace std;

// a^n mod を計算する
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
    long long m, n; cin >> m >> n;
    cout << modpow(m, n, 1000000007) << endl;
}
