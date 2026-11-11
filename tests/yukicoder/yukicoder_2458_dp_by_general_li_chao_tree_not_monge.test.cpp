#define PROBLEM "https://yukicoder.me/problems/no/2458"

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "general_li_chao_tree.hpp"

/*
    dp[j] = min_{j < i} (dp[i] - Q[i] × (Q[j]))
*/
int main() {
    const ll LIM = 210000, INF = 1LL<<60;
    ll N;
    cin >> N;
    vll Q(N);
    cin >> Q;

    auto make_func = [&](ll a, ll b) {
        return [&, a, b](ll q) -> ll {
            return a * q + b;
        };
    };
    DynamicGeneralLiChaoTree<ll, decltype(make_func(0, 0))> cht(-LIM, LIM);
    vll dp(N, INF);
    dp[0] = 0;
    cht.insert(make_func(-Q[0], dp[0]));
    REP2(j, 1, N) {
        dp[j] = cht.query(Q[j]).first;
        cht.insert(make_func(-Q[j], minll(dp[j], 0)));
    }
    ll res = -INF;
    REP2(i, 1, N) chmax(res, -dp[i]);
    cout << res << endl;
}
