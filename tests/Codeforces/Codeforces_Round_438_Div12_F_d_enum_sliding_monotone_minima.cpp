#define PROBLEM "https://codeforces.com/contest/868/problem/F"

#include <bits/stdc++.h>
using namespace std;

#include "monge_enum_d_shortest_paths.hpp"

// Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem
// by Monge 単一始点 d-辺最短路の d = 1, 2, ..., D における列挙 (by SMAWK 法, in O(ND))
int main() {
    long long N, K;
    cin >> N >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];

    vector<long long> nums(N + 1, 0);
    long long score = 0;
    auto add = [&](int i) -> void {
        score += nums[A[i]];
        nums[A[i]]++;
    };
    auto del = [&](int i) -> void {
        nums[A[i]]--;
        score -= nums[A[i]];
    };
    auto get = [&]() -> long long {
        return score;
    };
    MongeShortestPathWithDEdges<long long> msp;
    auto res = msp.solve(N, add, del, get, K, "monotone");  // we should use slide access
    cout << res[N][K].first << endl;
}
