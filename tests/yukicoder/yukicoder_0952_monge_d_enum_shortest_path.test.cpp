#define PROBLEM "https://yukicoder.me/problems/no/952"

#include <bits/stdc++.h>
using namespace std;

#include "monge_enum_d_shortest_paths.hpp"

// yukicoder No.952 危険な火薬庫
/*
    末尾にドアを追加する。
    ノード 0, 1, ..., N, N+1
    f(i, j) = (S[j-1] - S[i])^2
    選ぶドアが k 個 → 区間の個数は、N+1-k
*/
int main() {
    int N;
    cin >> N;
    vector<long long> A(N + 1, 0), S(N + 2, 0);
    for (int i = 0; i < N; i++) cin >> A[i], S[i + 1] = S[i] + A[i];
    auto f = [&](int i, int j) -> long long { 
        return (S[j - 1] - S[i]) * (S[j - 1] - S[i]);
    };
    MongeShortestPathWithDEdges<long long> msp;
    auto res = msp.solve(N + 1, f, N + 1);
    for (int k = 1; k <= N; k++) cout << res[N + 1][N + 1 - k].first << endl;
}
