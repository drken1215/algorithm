#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// 文字列 S, T の LCS を求める
string find_lcs(const string &S, const string &T) {
    vector<vector<int>> dp(S.size()+1, vector<int>(T.size()+1, 0));
    for (int i = 0; i < S.size(); ++i) {
        for (int j = 0; j < T.size(); ++j) {
            dp[i+1][j+1] = max(dp[i+1][j+1], dp[i][j] + (S[i] == T[j]));
            dp[i+1][j+1] = max(dp[i+1][j+1], dp[i+1][j]);
            dp[i+1][j+1] = max(dp[i+1][j+1], dp[i][j+1]);
        }
    }
    
    // recoustruct
    string res = "";
    int i = (int)S.size(), j = (int)T.size();
    while (i > 0 && j > 0) {
        if (dp[i][j] == dp[i-1][j]) --i;
        else if (dp[i][j] == dp[i][j-1]) --j;
        else {
            --i, --j;
            res = S[i] + res;
        }
    }
    return res;
}
