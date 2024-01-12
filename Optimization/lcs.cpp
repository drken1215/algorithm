//
// 最長共通部分列 (LCS) を求める
//
// verifed
//   EDPC F - LCS
//     https://atcoder.jp/contests/dp/tasks/dp_f
//


#include <bits/stdc++.h>
using namespace std;


template<class T> bool chmax(T& a, T b) { if (a < b) { a = b; return true; } return false; }
template<class T> bool chmin(T& a, T b) { if (a > b) { a = b; return true; } return false; }

string findLCS(const string &S, const string &T) {
    // DP
    vector<vector<int>> dp(S.size()+1, vector<int>(T.size()+1, 0));
    for (int i = 0; i < S.size(); ++i) {
        for (int j = 0; j < T.size(); ++j) {
            chmax(dp[i+1][j+1], dp[i][j] + (S[i] == T[j]));
            chmax(dp[i+1][j+1], dp[i+1][j]);
            chmax(dp[i+1][j+1], dp[i][j+1]);
        }
    }
    
    // recoustruct
    string res = "";
    int i = (int)S.size(), j = (int)T.size();
    while (i > 0 && j > 0)
    {
        if (dp[i][j] == dp[i-1][j]) --i;
        else if (dp[i][j] == dp[i][j-1]) --j; // DP の遷移を遡る
        else {
            --i, --j;
            res = S[i] + res;
        }
    }
    return res;
}



//------------------------------//
// Examples
//------------------------------//

void EDPC_F() {
    string S, T;
    cin >> S >> T;
    cout << findLCS(S, T) << endl;
}


int main() {
    EDPC_F();
}
