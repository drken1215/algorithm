#define PROBLEM "https://mofecoder.com/contests/monoxercon_202508/tasks/monoxercon_202508_k"

#include <bits/stdc++.h>
using namespace std;

#include "sparse_table.hpp"
#include "monge_enum_d_shortest_paths.hpp"

// Monoxer Programming Contest for Engineers K - Coupon
/*
    f(i, j) = (S[j] - S[i]) - max_{i ≦ k < j} A[k] / 2 + (j - i - 1)^2
*/
int main() {
    long long N;
    cin >> N;
    vector<long long> A(N), S(N+1, 0);
    for (int i = 0; i < N; i++) cin >> A[i], S[i+1] = S[i] + A[i];
    SparseTable<long long> sparsetable(A, [](long long p, long long q){return max(p, q);});
    auto f = [&](int i, int j) -> long long {
        return S[j] - S[i] - sparsetable.get(i, j)/2 + (j - i - 1) * (j - i - 1);
    };
    MongeShortestPathWithDEdges<long long> msp;
    auto res = msp.solve(N, f, N);
    for (int i = 1; i <= N; i++) cout << res[i] << endl;
}
