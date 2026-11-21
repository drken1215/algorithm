#define PROBLEM "https://judge.yosupo.jp/problem/min_plus_convolution_convex_convex"

#include <bits/stdc++.h>
using namespace std;

#include "min_plus_convolution_convex_convex.hpp"

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    long long N, M;
    cin >> N >> M;
    vector<long long> A(N), B(M);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < M; i++) cin >> B[i];
    auto res = MinPlusConvolutionConvexConvex(A, B);
    for (int i = 0; i < (int)res.size(); i++) cout << res[i] << " ";
    cout << endl;
}
