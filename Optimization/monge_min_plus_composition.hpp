#pragma once
#include <bits/stdc++.h>
#include <cassert>
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
