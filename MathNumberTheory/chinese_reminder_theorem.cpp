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
//   ACL Contest 1 B - Sum is Multiple
//     https://atcoder.jp/contests/acl1/tasks/acl1_b
//


/*
    x ≡ r[i] (mod. m[i]) の解を x ≡ R (mod. M) として {R, M} を求める
    解なしのときは {0, -1}
*/


#include <iostream>
#include <vector>
using namespace std;


long long extGcd(long long a, long long b, long long &p, long long &q) {
    if (b == 0) { 
        p = 1, q = 0; 
        return a; 
    }
    long long d = extGcd(b, a % b, q, p);
    q -= a / b * p;
    return d;
}

pair<long long, long long> ChineseRem(const vector<long long> &vr, const vector<long long> &vm) {
    if (vr.empty() || vm.empty()) return make_pair(0, 1);
    long long R = vr[0], M = vm[0];
    for (int i = 1; i < (int)vr.size(); ++i) {
        long long p, q, r = vr[i], m = vm[i];
        if (M < m) swap(M, m), swap(R, r); // prevent overflow
        long long d = extGcd(M, m, p, q); // p is inv of M/d (mod. m/d)
        if ((r - R) % d != 0) return make_pair(0, -1);
        long long md = m / d;
        long long tmp = (r - R) / d % md * p % md;
        R += M * tmp, M *= md;
    }
    R %= M;
    if (R < 0) R += M;
    return make_pair(R, M);
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
