#define PROBLEM "https://yukicoder.me/problems/no/952"

#include <bits/stdc++.h>
using namespace std;

#include "int128.hpp"
#include "monotone_convex_hull_trick.hpp"

// yukicoder No.952 危険な火薬庫
/*
    末尾にドアを追加する。
    ノード 0, 1, ..., N, N+1
    f(i, j) = (S[j-1] - S[i])^2
    選ぶドアが k 個 → 区間の個数は、N+1-k
*/
int main() {
    int N;
    cin >> N;
    N++;
    vector<i128> A(N, 0), S(N+1, 0);
    for (int i = 0; i < N; i++) {
        if (i) cin >> A[i];
        S[i+1] = S[i] + A[i];
    }

    vector chts(N+1, MonotoneConvexHullTrick<i128>());
    
    i128 INF = i128(1)<<100;
    vector dp(N+1, vector(N+1, INF));

    auto push = [&](int i, int j) -> void {
        if (i-1 >= 0) chts[i-j].insert(-S[i]*2, S[i]*S[i] + dp[i-1][j]);
    };

    dp[0][0] = 0;
    push(0, 0);

    for (int i = 1; i <= N; i++) {
        for (int j = 0; j < i; j++) {
            dp[i][j] = min(dp[i][j], dp[i-1][j]);
            dp[i][j] = min(dp[i][j], chts[i-j].query_monotone(S[i]) + S[i] * S[i]);
            push(i, j);
        }
    }

    for (int k = 1; k < N; k++) {
        cout << dp[N][k] << endl;
    }
}
