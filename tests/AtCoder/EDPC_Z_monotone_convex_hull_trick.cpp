#define PROBLEM "https://atcoder.jp/contests/dp/tasks/dp_z"

#include "monotone_convex_hull_trick.hpp"

// AtCoder EDPC Z - Frog 3
/*
    H は単調増加数列
    chmin(dp[j], dp[i] + (H[j] - H[i])^2 + C)
    i -> j のコスト：(H[j] - H[i])^2 ...... 差の凸関数は Monge
    スタート: 0, ゴール: N-1
*/
int main() {
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N);
    for (int i = 0; i < N; i++) cin >> H[i];

    const long long INF = 1LL<<60;
    const long long MAX = 1LL<<40;
    vector<long long> dp2(N, INF);
    /*
    　　dp[i] = min_j(dp[j] + (H[j] - H[i])² + C)
    　　dp2[i] = dp[i] + H[i]² とすると
    　　dp2[i] = min_j(-2 H[j] × H[i] + dp2[j]) + 2 H[i]² + C
    */
    dp2[0] = H[0] * H[0];
    CHT<long long> cht(INF);
    cht.insert(-H[0] * 2, dp2[0]);
    for (int i = 1; i < N; i++) {
        long long val = cht.query_monotone(H[i]);
        dp2[i] = min(dp2[i], val + H[i] * H[i] * 2 + C);
        cht.insert(-H[i] * 2, dp2[i]);
    }
    long long res = dp2[N-1] - H[N-1] * H[N-1];
    cout << res << endl;
}
