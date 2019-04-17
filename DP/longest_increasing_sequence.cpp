//
// LIS (Longest Increasing Sequence)
//   数列 a の最長増加部分列を求める
//     is_strong = true のとき狭義単調増加なもの、false のとき広義単調増加なもの
//
// verified
//   AOJ Course DPL_1_D Combinatorial - Longest Increasing Subsequence
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DPL_1_D&lang=jp
//


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


// dp[i] := 長さが i の増加部分列として最後尾の要素のとりうる最小値
template<class T> int LIS(vector<T> a,  bool is_strong = true) {
    const T INF = 1<<30; // to be set appropriately
    int n = (int)a.size();
    vector<T> dp(n, INF);
    for (int i = 0; i < n; ++i) {
        if (is_strong) *lower_bound(dp.begin(), dp.end(), a[i]) = a[i];
        else *upper_bound(dp.begin(), dp.end(), a[i]) = a[i];
    }
    return lower_bound(dp.begin(), dp.end(), INF) - dp.begin();
}


int main() {
    int N; cin >> N;
    vector<int> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    cout << LIS(a) << endl;
}
