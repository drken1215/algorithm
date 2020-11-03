//
// 離散対数
//
// cf.
//  「1000000007 で割ったあまり」の求め方を総特集！ 〜 逆元から離散対数まで 〜
//     https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a
//
// verified
//   m <= 500 まで単純解法と比較
//
//   ARC 042 D - あまり
//     https://atcoder.jp/contests/arc042/tasks/arc042_d
//


#include <bits/stdc++.h>
using namespace std;


// a^b
long long modpow(long long a, long long n, long long mod) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}

// a^-1
long long modinv(long long a, long long m) {
    long long b = m, u = 1, v = 0;
    while (b) {
        long long t = a / b;
        a -= t * b; swap(a, b);
        u -= t * v; swap(u, v);
    }
    u %= m;
    if (u < 0) u += m;
    return u;
}

// a^x ≡ b (mod. m) となる最小の正の整数 x を求める
long long modlog(long long a, long long b, int m) {
    a %= m, b %= m;

    // calc sqrt{M}
    long long lo = -1, hi = m;
    while (hi - lo > 1) {
        long long mid = (lo + hi) / 2;
        if (mid * mid >= m) hi = mid;
        else lo = mid;
    }
    long long sqrtM = hi;

    // {a^0, a^1, a^2, ..., a^sqrt(m)} 
    map<long long, long long> apow;
    long long amari = a;
    for (long long r = 1; r < sqrtM; ++r) {
        if (!apow.count(amari)) apow[amari] = r;
        (amari *= a) %= m;
    }

    // check each A^p
    long long A = modpow(modinv(a, m), sqrtM, m);
    amari = b;
    for (long long q = 0; q < sqrtM; ++q) {
        if (amari == 1 && q > 0) return q * sqrtM;
        else if (apow.count(amari)) return q * sqrtM + apow[amari];
        (amari *= A) %= m;
    }

    // no solutions
    return -1;
}



///////////////////////////////////////
// solver
///////////////////////////////////////

#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, vector<vector<T> > P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, set<T> P)
{ for(auto it : P) { s << "<" << it << "> "; } return s << endl; }
template<class T1, class T2> ostream& operator << (ostream &s, map<T1,T2> P)
{ for(auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s << endl; }

long long GCD(long long x, long long y) {
    if (y == 0) return x;
    else return GCD(y, x % y);
}

void subcheck(long long a, long long m) {
    vector<long long> correct(m, -1);
    long long p = a;
    for (int i = 1; i < m; ++i) {
        if (correct[p] == -1) correct[p] = i;
        p = (p * a) % m;
    }
    for (long long b = 1; b < m; ++b) {
        long long res = modlog(a, b, m);
        if (res != correct[b]) {
            cout << "error: case of " << a << ", " << b << ", " << m << endl;
        }
    }
}

void allcheck() {
    for (long long m = 2; m <= 1000; ++m) {
        if (m % 100 == 0) COUT(m);
        for (long long a = 1; a < m; ++a) {
            if (GCD(a, m) > 1) continue;
            subcheck(a, m);
        }
    }
}

long long ARC042D() {
    long long X, P, A, B;
    cin >> X >> P >> A >> B;
    long long r = modlog(X, 1, P);
    if (A == 0) return 1;
    if (B/r - (A-1)/r >= 1) return 1;

    A %= r, B %= r;
    if (B - A + 1 <= 1000000) {
        long long val = modpow(X, A, P);
        long long res = P;
        for (long long i = A; i <= B; ++i) {
            res = min(res, val);
            val = (val * X) % P;
        }
        return res;
    }
    else {
        for (long long b = 1; b < P; ++b) {
            long long exp = modlog(X, b, P);
            if (A <= exp && exp <= B) return b;
        }
    }
    return P;
}

int main() {
    /*
    COUT(modlog(4, 1, 5));
    COUT(modlog(10, 10, 11));
    COUT(modlog(2, 1, 7));
    allcheck();
    */
    
    cout << ARC042D() << endl;
}
