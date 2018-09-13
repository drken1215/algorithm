//
// 素因数分解
//
// verified
//   ARC 017 A - 素数、コンテスト、素数
//     https://beta.atcoder.jp/contests/arc017/tasks/arc017_1
//

/*
    n が素数かどうかを O(√n) で判定, n は 0 以上の整数とする
*/


#include <iostream>
using namespace std;

bool is_prime(long long n) {
    if (n <= 1) return false;
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p == 0) return false;
    }
    return true;
}

int main() {
    int n; cin >> n;
    if (is_prime(n)) cout << "YES" << endl;
    else cout << "NO" << endl;
}
