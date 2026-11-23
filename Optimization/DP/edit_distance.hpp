#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// 文字列 S, T の編集距離を求める
int find_edit_distance(const string &S, const string &T) {
    vector<vector<int>> dp(S.size()+1, vector<int>(T.size()+1, S.size()+T.size()));
    dp[0][0] = 0;
    for (int i = 0; i < S.size(); ++i) dp[i+1][0] = i+1;
    for (int j = 0; j < T.size(); ++j) dp[0][j+1] = j+1;
    for (int i = 0; i < S.size(); ++i) {
        for (int j = 0; j < T.size(); ++j) {
            dp[i+1][j+1] = min(dp[i+1][j+1], dp[i][j] + (S[i] != T[j]));
            dp[i+1][j+1] = min(dp[i+1][j+1], dp[i+1][j] + 1);
            dp[i+1][j+1] = min(dp[i+1][j+1], dp[i][j+1] + 1);
        }
    }
    return dp[S.size()][T.size()];
}
