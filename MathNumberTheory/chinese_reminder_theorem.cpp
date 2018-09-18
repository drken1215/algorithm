//
// 中国剰余定理
//
// cf.
//   中国剰余定理 (CRT) の解説と、それを用いる問題のまとめ
//     https://qiita.com/drken/items/ae02240cd1f8edfc86fd
//
// verified
//   AOJ 2659 箸
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2659
//


/*
    x ≡ b[i] (mod. m[i]) の解を x ≡ B (mod. M) として {B, M} を求める
    解なしのときは {0, -1}
*/


#include <iostream>
#include <vector>
using namespace std;

inline long long mod(long long a, long long m) {
    return (a % m + m) % m;
}

long long extGcd(long long a, long long b, long long &p, long long &q) {
    if (b == 0) { p = 1; q = 0; return a; }
    long long d = extGcd(b, a%b, q, p);
    q -= a/b * p;
    return d;
}

pair<long long, long long> ChineseRem(const vector<long long> &b, const vector<long long> &m) {
    long long r = 0, M = 1;
    for (int i = 0; i < (int)b.size(); ++i) {
        long long p, q;
        long long d = extGcd(M, m[i], p, q); // p is inv of M/d (mod. m[i]/d)
        if ((b[i] - r) % d != 0) return make_pair(0, -1);
        long long tmp = (b[i] - r) / d * p % (m[i]/d);
        r += M * tmp;
        M *= m[i]/d;
    }
    return make_pair(mod(r, M), M);
}



int main() {
    long long N;
    int M, D;
    cin >> N >> M >> D;
    vector<long long> A(M);
    for (int i = 0; i < M; ++i) cin >> A[i];
    bool ok = true;
    for (int i = 0; i < D; ++i) {
        vector<long long> b, m;
        for (int j = 0; j < M; ++j) {
            int R; cin >> R;
            if (R != -1) b.push_back(R), m.push_back(A[j]);
        }
        if (b.empty()) continue;
        
        pair<long long, long long> tmp = ChineseRem(b, m);
        
        if (tmp.second == -1) ok = false;
        if (N < tmp.first) ok = false;
        
        N = N - (N - tmp.first) % tmp.second;
    }
    if (ok) cout << N << endl;
    else cout << -1 << endl;
}
