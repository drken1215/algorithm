//
// ヒストグラムにおいて、左右にどこまで伸ばせるかを求める
//
// verifed
//   AOJ Course DPL_3_C - ヒストグラムの中の最大長方形
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DPL_3_C&lang=ja
//
//   AOJ Couser DPL_3_B - 最大長方形
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DPL_3_B&lang=ja
//
//   ABC 311 G - One More Grid Task
//     https://atcoder.jp/contests/abc311/tasks/abc311_g
//


#include <bits/stdc++.h>
using namespace std;


// min value of H within [left[i], right[i]) <= H[i]
pair<vector<int>, vector<int>> solve_left_right(const vector<long long> &H) {
    int N = (int)H.size();
    vector<int> left(N, 0), right(N, N);
    
    // left
    stack<pair<long long, int>> stack_left;
    for(int i = 0; i < N; ++i) {
        while (!stack_left.empty() && H[i] <= stack_left.top().first)
            stack_left.pop();
        if (!stack_left.empty()) left[i] = stack_left.top().second + 1;
        stack_left.push({H[i], i});
    }
    
    // right
    stack<pair<long long, int>> stack_right;
    for(int i = N-1; i >= 0; --i) {
        while (!stack_right.empty() && H[i] <= stack_right.top().first)
            stack_right.pop();
        if (!stack_right.empty()) right[i] = stack_right.top().second;
        stack_right.push({H[i], i});
    }
    return {left, right};
}

// max area of rectangle in histogram
long long max_area_in_histogram(const vector<long long> &H) {
    auto [left, right] = solve_left_right(H);
    long long res = 0;
    for (int i = 0; i < H.size(); ++i) {
        res = max(res, H[i] * (right[i] - left[i]));
    }
    return res;
}



//------------------------------//
// Examples
//------------------------------//

void AOJ_Course_DPL_3_C() {
    int N;
    cin >> N;
    vector<long long> H(N);
    for (int i = 0; i < N; ++i) cin >> H[i];
    cout << max_area_in_histogram(H) << endl;
}

void AOJ_Course_DPL_3_B() {
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

void ABC_311_G() {
    // 入力
    int H, W;
    cin >> H >> W;
    vector<vector<long long>> A(H, vector<long long>(W));
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) cin >> A[i][j];
    
    // 2 次元累積和
    vector<vector<long long>> S(H+1, vector<long long>(W+1, 0));
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            S[i+1][j+1] = S[i+1][j] + S[i][j+1] - S[i][j] + A[i][j];
        }
    }
    
    // 長方形領域 [lx, rx) x [ly, ry) の総和
    auto calc = [&](int lx, int rx, int ly, int ry) -> long long {
        return S[min(rx, H)][min(ry, W)]
        - S[min(rx, H)][ly] - S[lx][min(ry, W)] + S[lx][ly];
    };
    
    // 最小値が v 以上の場合を求める (v 未満のマスを禁止する)
    long long res = 0;
    for (int v = 0; v <= 300; ++v) {
        // 各マスを起点として下側に可能マスが何個連続するかを求める
        vector<vector<long long>> len(H+1, vector<long long>(W, 0));
        for (int j = 0; j < W; ++j) {
            for (int i = H-1; i >= 0; --i) {
                if (A[i][j] < v) len[i][j] = 0;
                else len[i][j] = len[i+1][j] + 1;
            }
        }
        
        // 行 i を起点として、その下側のヒストグラムを考える
        for (int i = 0; i < H; ++i) {
            auto [left, right] = solve_left_right(len[i]);
            
            // マス (i, j) を起点として、下方向、左方向、右方向に限界まで伸ばした長方形領域の総和
            for (int j = 0; j < W; ++j) {
                long long sum = calc(i, i+len[i][j], left[j], right[j]);
                res = max(res, sum * v);
            }
        }
    }
    cout << res << endl;
}


int main() {
    //AOJ_Course_DPL_3_C();
    //AOJ_Course_DPL_3_B();
    ABC_311_G();
}

