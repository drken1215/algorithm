//
// 白黒グリッド内の、白色のみからなる最大正方形を求める
//
// verifed
//   AOJ Course DPL_3_A - 最大正方形
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DPL_3_A&lang=ja
//
//   AtCoder ABC 311 E - Defect-free Squares
//     https://atcoder.jp/contests/abc311/tasks/abc311_e
//


#include <bits/stdc++.h>
using namespace std;


void AOJ_DPL_3_A() {
    // 0 is ok, 1 is ng
    int H, W;
    cin >> H >> W;
    vector<vector<int>> C(H, vector<int>(W));
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) cin >> C[i][j];
    
    int res = 0;
    
    // dp[i+1][j+1] := マス (i, j) を右下マスとする最大の正方形の一辺の長さ
    vector<vector<int>> dp(H+1, vector<int>(W+1, 0));
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            if (C[i][j] == 0) {
                dp[i+1][j+1] = min({dp[i+1][j], dp[i][j+1], dp[i][j]}) + 1;
                res = max(res, dp[i+1][j+1]);
            }
        }
    }
    cout << res*res << endl;
}

void ABC_311_E() {
    // 0 is ok, 1 is ng
    int H, W, N;
    cin >> H >> W >> N;
    vector<vector<int>> S(H, vector<int>(W, 0));
    for (int i = 0; i < N; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        S[a][b] = 1;
    }
    
    // dp[i+1][j+1] := マス (i, j) を右下マスとする最大の正方形の一辺の長さ
    vector<vector<int>> dp(H+1, vector<int>(W+1, 0));
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            if (S[i][j] == 0) {
                dp[i+1][j+1] = min({dp[i+1][j], dp[i][j+1], dp[i][j]}) + 1;
            }
        }
    }
    
    // 集計
    long long res = 0;
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) res += dp[i+1][j+1];
    cout << res << endl;
}


int main() {
    //AOJ_DPL_3_A();
    ABC_311_E();
}

