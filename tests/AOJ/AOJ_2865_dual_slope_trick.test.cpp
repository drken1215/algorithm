#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/problems/2865"

#include <bits/stdc++.h>
using namespace std;

#include "dual_slope_trick.hpp"

// AOJ 2865 Farm Village
/*
    dp[x] := 左から i 個分のノードから溢れているフロー量が x であるときの最小コスト (x の分も含む)

    nex[x] = min(dp[x+1], dp[x] + G[i], dp[x-1] + 2G[i]) + D[i]|x|
    g(y) = G[i]y (-1 ≦ y ≦ 1), ∞ (それ以外) として、

    nex[x] = min_{y}(dp[x - y] + g(y)) + D[i]|x| + G[i]
*/
int main() {
    long long N;
    cin >> N;
    vector<long long > D(N-1), G(N);
    for (int i = 0; i < N-1; i++) cin >> D[i];
    for (int i = 0; i < N; i++) cin >> G[i];
    DualSlopeTrick<long long, long long> dp(1LL << 45);
    for (int i = 0; i < N; i++) {
        DualSlopeTrick<long long, long long> g(1LL << 45);
        g.pushL(G[i]), g.pushR(G[i]);
        dp.min_plus_convolution(g);
        dp.add_abs(D[i], 0);
        dp.add_const(G[i]);
    }
    long long res = dp.get_f0();
    cout << res << endl;
}
