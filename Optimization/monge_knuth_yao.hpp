#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// f: Monge, interval monotone ((N+1) x (N+1))
// dp[i][j] = min_{i <= k < j} (dp[i][k] + dp[k][j]) + f(i, j), in O(N^2)
template<class VAL, class FUNC> VAL KnuthYao(int N, const FUNC &f) {
    VAL INF = numeric_limits<VAL>::max() / 2;
    vector<vector<VAL>> dp(N+1, vector<VAL>(N+1, INF));
    vector<vector<int>> K(N+1, vector<int>(N+1, -1));

    for (int i = 0; i <= N; i++) {
        dp[i][i] = 0, K[i][i] = i;
        if (i < N) dp[i][i+1] = f(i, i+1), K[i][i+1] = i;
    }
    for (int between = 2; between <= N; between++) {
        for (int i = 0; i + between <= N; i++) {
            int j = i + between;
            for (int k = K[i][j-1]; k <= K[i+1][j]; k++) {
                VAL tmp = dp[i][k] + dp[k][j] + f(i, j);
                if (dp[i][j] > tmp) {
                    dp[i][j] = tmp;
                    K[i][j] = k;
                }   
            }
        }
    }
    return dp[0][N];
}
