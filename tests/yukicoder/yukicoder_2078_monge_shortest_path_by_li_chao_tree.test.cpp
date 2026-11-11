#define PROBLEM "https://yukicoder.me/problems/no/2078"

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "monge_shortest_path_by_cht.hpp"

/*
    dp[j] = min_{j < i} (dp[i] + func(i, j))
    func(i, j) = (j - i)^2 + X
*/
void solve() {
    ll X, A;
    cin >> X >> A;
    auto func = [&](ll i, ll j) -> ll { return (j - i) * (j - i) + X; };
    auto res = MongeShortestPathByCHT<ll>(A, func);
    cout << res.back().first << endl;
}

int main() {
    int T;
    cin >> T;
    while (T--) solve();
}
