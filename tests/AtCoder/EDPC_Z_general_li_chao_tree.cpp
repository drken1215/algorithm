#define PROBLEM "https://atcoder.jp/contests/dp/tasks/dp_z"

#include <bits/stdc++.h>
using namespace std;

#include "macro.hpp"
#include "general_li_chao_tree.hpp"

// AtCoder EDPC Z - Frog 3
/*
    H は単調増加数列
    chmin(dp[j], dp[i] + (H[j] - H[i])^2 + C)
    i -> j のコスト：(H[j] - H[i])^2 ...... 差の凸関数は, i, j 全域で Monge
    スタート: 0, ゴール: N-1
*/
int main() {
    const long long INF = LL_INF;
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N), dp(N, INF);
    for (long long i = 0; i < N; i++) cin >> H[i];
    auto func = [&](int i, int j) -> long long {
        return (H[j] - H[i]) * (H[j] - H[i]) + C;
    };
    auto make_func = [&](int i) {
        return [&, i](long long j) -> long long {
            return dp[i] + func(i, j);
        };
    };
    dp[0] = 0;
    GeneralLiChaoTree<long long, decltype(make_func(0))> cht(1, N);
    for (int j = 1; j < N; j++) {
        cht.insert(make_func(j - 1));
        dp[j] = cht.query(j).first;
    }
    cout << dp[N - 1] << endl;
}
