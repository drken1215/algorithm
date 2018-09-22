//
// 文字列 DP で頻繁に必要になる前処理
//
// verified:
//   ARC 081 E - Don't Be a Subsequence
//     https://beta.atcoder.jp/contests/arc081/tasks/arc081_c
//

/*
    res[i][c] := i 文字目以降で最初に文字 c が登場する index (存在しないときは n)
 */


#include <iostream>
#include <string>
#include <vector>
using namespace std;


// res[i][c] := i 文字目以降で最初に文字 c が登場する index (存在しないときは n)
vector<vector<int> > calcNext(const string &S) {
    int n = (int)S.size();
    vector<vector<int> > res(n+1, vector<int>(26, n));
    for (int i = n-1; i >= 0; --i) {
        for (int j = 0; j < 26; ++j) res[i][j] = res[i+1][j];
        res[i][S[i]-'a'] = i;
    }
    return res;
}



// chmin
template<class T> inline bool chmin(T& a, T b) { if (a > b) { a = b; return 1; } return 0; }

int main() {
    string S; cin >> S;
    int n = (int)S.size();
    
    // 前処理配列
    auto next = calcNext(S);
    
    // DP
    vector<int> dp(n+1, 1<<29); // DP テーブル
    vector<pair<char, int> > recon(n+1, {'?', n}); // 復元用テーブル
    dp[n] = 1;
    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j < 26; ++j) {
            // 次の文字がないとき
            if (next[i][j] == n) {
                if (dp[i] > 1)
                    dp[i] = 1, recon[i] = {'a'+j, n};
            }
            // 次の文字があるとき
            else if (chmin(dp[i], dp[next[i][j] + 1] + 1))
                recon[i] = {'a'+j, next[i][j] + 1};
        }
    }
    
    // 復元
    string res = "";
    int index = 0;
    while (index < n) {
        auto p = recon[index];
        res += p.first;
        index = p.second;
    }
    cout << res << endl;
}
