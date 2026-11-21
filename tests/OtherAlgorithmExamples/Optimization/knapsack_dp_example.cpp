//
// ナップサック問題を解く DP
//
// reference:
//   drken (Qiita): 典型的な DP (動的計画法) のパターンを整理 Part 1 ～ ナップサック DP 編 ～
//     https://qiita.com/drken/items/a5e6fe22863b7992efdb
//
// verified:
//   EDPC D - Knapsack 1
//     https://atcoder.jp/contests/dp/tasks/dp_d
//


#include <iostream>
#include <vector>
using namespace std;

template<class T> void chmax(T &a, T b) {
    if (a < b) a = b;
}

int main() {
    int N, W;
    cin >> N >> W;
    vector<int> weight(N), value(N);
    for (int i = 0; i < N; ++i) {
        cin >> weight[i] >> value[i];
    }

    // dp 配列
    vector<long long> dp(W + 1, 0);

    // ループ
    for (int i = 0; i < N; ++i) {
        // 更新後の配列を用意する
        vector<long long> nex(W + 1, 0);

        // DP の 1 ステップ
        for (int w = 0; w <= W; ++w) {
            chmax(nex[w], dp[w]);
            if (w + weight[i] <= W) {
                chmax(nex[w + weight[i]], dp[w] + value[i]);
            }
        }

        // dp と nex を swap する
        swap(dp, nex);
    }

    // 答え
    cout << dp[W] << endl;
}
