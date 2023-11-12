//
// 編集距離を求める
//
// verifed
//   典型アルゴリズム問題集 上級〜エキスパート編 B - 編集距離
//     https://atcoder.jp/contests/pastbook2022/tasks/pastbook2022_b
//


#include <bits/stdc++.h>
using namespace std;

template<class T> bool chmax(T& a, T b) { if (a < b) { a = b; return true; } return false; }
template<class T> bool chmin(T& a, T b) { if (a > b) { a = b; return true; } return false; }

int EditDistance(const string &S, const string &T) {
    vector<vector<int>> dp(S.size()+1, vector<int>(T.size()+1, S.size()+T.size()));
    dp[0][0] = 0;
    for (int i = 0; i < S.size(); ++i) dp[i+1][0] = i+1;
    for (int j = 0; j < T.size(); ++j) dp[0][j+1] = j+1;
    for (int i = 0; i < S.size(); ++i) {
        for (int j = 0; j < T.size(); ++j) {
            chmin(dp[i+1][j+1], dp[i][j] + (S[i] != T[j]));
            chmin(dp[i+1][j+1], dp[i+1][j] + 1);
            chmin(dp[i+1][j+1], dp[i][j+1] + 1);
        }
    }
    return dp[S.size()][T.size()];
}



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void PAST_Hard_Expert_B() {
    int M, N;
    string S, T;
    cin >> M >> N >> S >> T;
    cout << EditDistance(S, T) << endl;
}


int main() {
    PAST_Hard_Expert_B();
}

