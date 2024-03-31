//
// m で割って r 余る, x 以上の最小の整数
//
// verified:
//   ABC 129 F - Takahashi's Basics in Education and Learning
//     https://atcoder.jp/contests/abc129/tasks/abc129_f
//
//   AGC 066 A - Adjacent Difference
//     https://atcoder.jp/contests/agc066/tasks/agc066_a
//


#include <bits/stdc++.h>
using namespace std;

// m で割って r 余る x 以上の最小の整数
long long amari_lower_bound(long long x, long long m, long long r) {
    x -= r;
    long long rx = (x % m + m) % m;
    return (rx > 0 ? x + (m - rx) + r : x + r);
}



//------------------------------//
// Examples
//------------------------------//

void AGC_066_A() {
    int N, D;
    cin >> N >> D;
    vector A(N, vector(N, 0));
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) cin >> A[i][j];
    
    for (int k = 0; k < D * 2; ++k) {
        long long cost = 0;
        auto res = A;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                int r = ((i + j) % 2 == 0 ? k : (k + D) % (D * 2));
                long long up = amari_lower_bound(A[i][j], D * 2, r);
                long long down = up - D * 2;
                
                if (abs(up - A[i][j]) < abs(down - A[i][j])) {
                    cost += abs(up - A[i][j]);
                    res[i][j] = up;
                } else {
                    cost += abs(down - A[i][j]);
                    res[i][j] = down;
                }
            }
        }
        if (cost * 2 <= N * N * D) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    cout << res[i][j] << " ";
                }
                cout << endl;
            }
            return;
        }
    }
    
}


int main() {
    AGC_066_A();
}

