#define PROBLEM "https://atcoder.jp/contests/abc250/tasks/abc250_g"

#include <bits/stdc++.h>
using namespace std;

#include "dual_slope_trick.hpp"

// ABC 250 G - Stonks
/*
    dp[x] := x 個持っている状態でのお金の減少分の最小値（負になると嬉しい）

    nex[x] = min(dp[x], dp[x-1] + P[i], dp[x+1] - P[i]) (x >= 0)
    → nex[x] = min_{-1 ≦ y ≦ 1} (dp[x - y] + P[i]y) (x >= 0)
*/
int main() {
    long long N;
    cin >> N;
    vector<long long> P(N);
    for (int i = 0; i < N; i++) cin >> P[i];
    DualSlopeTrick<long long, long long> dp(1LL << 45);
    for (int i = 0; i < N; i++) {
        DualSlopeTrick<long long, long long> g(1LL << 45);
        g.pushL(P[i]), g.pushR(P[i]);
        dp.min_plus_convolution(g);
        dp.clear_left();
    }
    long long res = -dp.get_f0();
    cout << res << endl;
}
