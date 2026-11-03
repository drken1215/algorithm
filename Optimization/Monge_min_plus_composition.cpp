//
// 上三角 (N+1) x (N+1) Monge 行列 A, B に対して
//   C[i][j] = min_{i <= k < j} (A[i][k] + A[k][j])
//
// verified
//


#include <bits/stdc++.h>
using namespace std;


// C[i][j] = min_{i <= k < j} (A[i][k] + B[k][j])
// A, B: 上三角 (N+1) x (N+1) Monge
template<class VAL, class FUNCA, class FUNCB> vector<vector<VAL>> min_plus_composition
(int N, const FUNCA &A, const FUNCB &B) {
    VAL INF = numeric_limits<VAL>::max() / 2;
    vector<vector<VAL>> res(N+1, vector<VAL>(N+1, INF));
    vector<int> argmin(N+1);
    for (int i = 0; i <= N; i++) res[i][i] = A(i, i) + B(i, i), argmin[i] = i;
    for (int s = 1; s <= N; s++) {
        vector<int> argmin2(N + 1 - s);
        for (int i = 0; i <= N - s; i++) {
            int j = i + s, p = argmin[i], q = argmin[i + 1];
            for (int k = p; k <= q; k++) {
                if (res[i][j] > A(i, k) + B(k, j)) {
                    res[i][j] = A(i, k) + B(k, j);
                    argmin2[i] = k;
                }
            }
        }
        swap(argmin, argmin2);
    }
    return res;
}


//------------------------------//
// Examples
//------------------------------//

void mini_test() {
    int N = 5;
    auto A = [&](int i, int j) -> int { return (j - i) * (j - i); };
    auto B = [&](int i, int j) -> int { return max(j - i, i - j) * 3; };
    auto res = min_plus_composition<int>(N, A, B);
    for (int i = 0; i <= N; i++) {
        for (int j = i; j <= N; j++) cout << i << ", " << j << ": " << res[i][j] << endl;
    }
}


int main() {
    mini_test();
}