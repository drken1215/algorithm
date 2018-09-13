//
// 確率的素数判定 (Miller-Rabin)
//

/*
    ある程度巨大な整数 n に対してもそこそこ高い確率で正しく素数判定する。
    なお、n < 341,550,071,728,321については正しいことが保証される
*/


#include <iostream>
using namespace std;


inline long long mod(long long a, long long m) {
    return (a % m + m) % m;
}

long long mul(long long a, long long b, long long m) {
    if (b == 0) return 0;
    long long res = mul(mod(a + a, m), b>>1, m);
    if (b & 1) res = mod(res + a, m);
    return res;
}

long long pow(long long a, long long n, long long m) {
    if (n == 0) return 1 % m;
    long long t = pow(a, n/2, m);
    t = mul(t, t, m);
    if (n & 1) t = mul(t, a, m);
    return t;
}

bool is_prime(long long n) {
    if (n == 2) return true;
    if (!(n & 1) || n <= 1) return false;
    int a[] = {2, 3, 5, 7, 11, 13, 17, -1};
    long long s = 0, d = n-1;
    while (!(d & 1)) { ++s; d = d>>1; }
    for (int i = 0; a[i] != -1 && a[i] < n; ++i) {
        long long x = pow(a[i], d, n);
        int t;
        if (x != 1) {
            for (t = 0; t < s; ++t) {
                if (x == n-1) break;
                x = mul(x, x, n);
            }
            if (t == s) return false;
        }
    }
    return true;
}


int main() {
    int n; cin >> n;
    if (is_prime(n)) cout << "YES" << endl;
    else cout << "NO" << endl;
}
