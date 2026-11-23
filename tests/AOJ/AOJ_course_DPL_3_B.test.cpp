#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/problems/DPL_3_B"

#include <bits/stdc++.h>
using namespace std;

#include "histogram.hpp"

int main() {
    int H, W;
    cin >> H >> W;
    vector<vector<int>> C(H, vector<int>(W));
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) cin >> C[i][j];
    
    // 各マスを起点として下側に何個 0 が連続するかを求める
    vector<vector<long long>> len(H+1, vector<long long>(W, 0));
    for (int j = 0; j < W; ++j) {
        for (int i = H-1; i >= 0; --i) {
            if (C[i][j] == 1) len[i][j] = 0;
            else len[i][j] = len[i+1][j] + 1;
        }
    }
    
    // 各行を起点として、下側に伸びる 0 の長さをヒストグラムに見立てて、最大長方形を求める
    long long res = 0;
    for (int i = 0; i < H; ++i) {
        res = max(res, max_area_in_histogram(len[i]));
    }
    cout << res << endl;
}
