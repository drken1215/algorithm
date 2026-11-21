#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// Min-Plus Convolution (Convex and Convex)
template<class T> vector<T> MinPlusConvolutionConvexConvex(vector<T> &A, vector<T> &B) {
    int N = (int)A.size(), M = (int)B.size();
    if (N == 0 || M == 0) return {};
    vector<T> res(N + M - 1);
    int a = 0, b = 0;
    res[0] = A[0] + B[0];
    for (int i = 1; i < N + M - 1; i++) {
        if (b == M - 1 || (a < N - 1 && A[a + 1] + B[b] < A[a] + B[b + 1])) {
            res[i] = A[++a] + B[b];
        } else {
            res[i] = A[a] + B[++b];
        }
    }
    return res;
}
