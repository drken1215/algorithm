#define PROBLEM "https://atcoder.jp/contests/abc218/tasks/abc218_h"

#include <bits/stdc++.h>
using namespace std;

#include "alien_dp.hpp"

int main() {
    long long N, R;
    cin >> N >> R;
    vector<long long> A(N - 1);
    for (int i = 0; i < N - 1; i++) cin >> A[i];
    auto cost = [&](int i, int j) -> long long {
        if (j - i <= 1) return 0;
        if (i == 0) return -A[j - 2];
        else return -A[i - 1] - A[j - 2];
    };
    AlienDP<long long> ad;
    //ad.debug(N, cost);
    auto [resR, cntR] = ad.solve(N, cost, R);
    auto [resB, cntB] = ad.solve(N, cost, N - R);
    auto res = max(-resR, -resB);
    cout << res << endl;
}
