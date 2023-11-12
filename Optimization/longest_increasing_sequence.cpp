//
// LIS (Longest Increasing Sequence)
//   ���� a �κ�Ĺ������ʬ������
//     is_strong = true �ΤȤ�����ñĴ���äʤ�Ρ�false �ΤȤ�����ñĴ���äʤ��
//
// verified
//   AOJ Course DPL_1_D Combinatorial - Longest Increasing Subsequence
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DPL_1_D&lang=jp
//


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


// dp[i] := Ĺ���� i ��������ʬ��Ȥ��ƺǸ��������ǤΤȤꤦ��Ǿ���
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
