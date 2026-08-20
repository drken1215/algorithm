//
// 上三角 (N+1) x (N+1) Monge 行列 A, B に対して
//   C[i][j] = min_{i <= k < j} (A[i][k] + A[k][j])
//
// verified
//   ABC 254 B - Mongeness
//     https://atcoder.jp/contests/abc224/tasks/abc224_b
//


#include <bits/stdc++.h>
using namespace std;


// check whether N x M matrix f is Monge or not, in O(NM)
template<class FUNC> bool is_monge
(int N, int M, const FUNC &f, bool upper_triangular = false, bool output_detail = true) {
    assert(N > 0 && M > 0);
    for (int i = 0; i + 1 < N; i++) {
        for (int j = (upper_triangular ? i+2 : 0); j + 1 < M; j++) {
            auto f00 = f(i, j), f01 = f(i, j + 1), f10 = f(i + 1, j), f11 = f(i + 1, j + 1);
            if (f00 + f11 > f10 + f01) {
                if (output_detail) {
                    cout << "Not Monge! " << endl;
                    cout << "f(" << i << ", " << j << ") = " << f(i, j) << ", ";
                    cout << "f(" << i << ", " << j + 1 << ") = " << f(i, j + 1) << endl;
                    cout << "f(" << i + 1 << ", " << j << ") = " << f(i + 1, j) << ", ";
                    cout << "f(" << i + 1 << ", " << j + 1 << ") = " << f(i + 1, j + 1) << endl;
                }
                return false;
            }
        }
    }
    return true;
}

// check whether N x M matrix f is anti-Monge or not, in O(NM)
template<class FUNC> bool is_anti_monge
(int N, int M, const FUNC &f, bool upper_triangular = false, bool output_detail = true) {
    assert(N > 0 && M > 0);
    for (int i = 0; i + 1 < N; i++) {
        for (int j = (upper_triangular ? i+2 : 0); j + 1 < M; j++) {
            auto f00 = f(i, j), f01 = f(i, j + 1), f10 = f(i + 1, j), f11 = f(i + 1, j + 1);
            if (f00 + f11 < f10 + f01) {
                if (output_detail) {
                    cout << "Not anti-Monge! " << endl;
                    cout << "f(" << i << ", " << j << ") = " << f(i, j) << ", ";
                    cout << "f(" << i << ", " << j + 1 << ") = " << f(i, j + 1) << endl;
                    cout << "f(" << i + 1 << ", " << j << ") = " << f(i + 1, j) << ", ";
                    cout << "f(" << i + 1 << ", " << j + 1 << ") = " << f(i + 1, j + 1) << endl;
                }
                return false;
            }
        }
    }
    return true;
}


//------------------------------//
// Examples
//------------------------------//

// ABC 254 B - Mongeness
void ABC_254_B() {
    int H, W;
    cin >> H >> W;
    vector A(H, vector(W, 0LL));
    for (int i = 0; i < H; i++) for (int j = 0; j < W; j++) cin >> A[i][j];
    auto f = [&](int i, int j) -> long long { return A[i][j]; };
    if (is_monge(H, W, f, false, false)) cout << "Yes" << endl;
    else cout << "No" << endl;
}


int main() {
    ABC_254_B();
}