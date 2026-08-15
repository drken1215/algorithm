//
// 累積和
//
// cf.
//   累積和を何も考えずに書けるようにする！
//    https://qiita.com/drken/items/56a6b68edef8fc605821
//
// verified
//   AOJ 0516 最大の和
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=0516
//


#include <iostream>
#include <vector>
using namespace std;
const long long INF = 1LL<<60; // 仮想的な無限大の値

int main() {
    // 入力
    int N, K;
    while (cin >> N >> K) {
        if (N == 0) break;
        vector<long long> a(N);
        for (int i = 0; i < N; ++i) cin >> a[i];

        // 累積和を前計算
        vector<int> s(N+1, 0);
        for (int i = 0; i < N; ++i) s[i+1] = s[i] + a[i];

        // 答えを求める
        long long res = -INF; // 最初は無限小の値に初期化しておく
        for (int i = 0; i <= N-K; ++i) {
            long long val = s[K+i] - s[i];
            if (res < val) res = val;
        }
        cout << res << endl;
    }
}
