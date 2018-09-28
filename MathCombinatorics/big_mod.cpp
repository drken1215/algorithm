//
// 巨大な mod 対策 (普通にやるとオーバーフローしてしまう)
//
// cf.
//   【ライブラリ】mod の値が大きいときの mod 演算
//     http://drken1215.hatenablog.com/entry/2018/02/08/113249
//
// verified:
//   AOJ 2353 Four Arithmetic Operations
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2353
//


/*
    a を b 回足しあげる、ただし繰り返し二乗法を流用する)
    という方法をとる
 
    
    問題例
    ・AOJ 2353 Four Arithmetic Operations
    ・CS Academy 068 DIV2 E - Sliding Product Sum
 */


#include <iostream>
using namespace std;


inline long long mod(long long a, long long m) {
    return (a % m + m) % m;
}

inline long long mul(long long a, long long b, long long m) {
    a = mod(a, m); b = mod(b, m);
    if (b == 0) return 0;
    long long res = mul(mod(a + a, m), b>>1, m);
    if (b & 1) res = mod(res + a, m);
    return res;
}

inline long long inv(long long a, long long m) {
    long long b = m, u = 1, v = 0;
    while (b) {
        long long t = a/b;
        a = mod(a - mul(t, b, m), m); swap(a, b);
        u = mod(u - mul(t, v, m), m); swap(u, v);
    }
    return mod(u, m);
}

long long pow(long long a, long long n, long long m) {
    if (n == 0) return 1 % m;
    long long t = pow(a, n/2, m);
    t = mul(t, t, m);
    if (n & 1) t = mul(t, a, m);
    return t;
}

// 級数和
long long ser(long long a, long long n, long long m) {
    if (n == 0) return 0;
    long long res = ser(mul(a, a, m), n/2, m);
    res = mul(res, a+1, m);
    if (n & 1) res = mod(mul(res, a, m) + 1, m);
    return res;
}



int main() {
    const long long MOD = 67280421310721LL;
    long long res = 0;
    int N; cin >> N;
    for (int i = 0; i < N; ++i) {
        long long o, y; cin >> o >> y;
        if (o == 1) res = mod(res + y, MOD);
        else if (o == 2) res = mod(res - y, MOD);
        else if (o == 3) res = mul(res, y, MOD);
        else res = mul(res, inv(y, MOD), MOD);
    }
    if (res <= 1LL<<31) cout << res << endl;
    else cout << res - MOD << endl;
}
