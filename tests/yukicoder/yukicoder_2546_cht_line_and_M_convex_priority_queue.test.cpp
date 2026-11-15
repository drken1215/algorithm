#define PROBLEM "https://yukicoder.me/problems/no/2546"

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "monotone_convex_hull_trick.hpp"

/*
    f[a,d](m) = (A*2 + D(m-1)) * m / 2
              = (-Dm + d-2A) * (-m/2)
*/
int main() {
    ll N, M, a, d, INF = 1LL << 60;
    cin >> N >> M;
    vll PA, PD;
    priority_queue<pll> que;
    REP(i, N) {
        cin >> a >> d;
        if (d >= 0) PA.EB(a), PD.EB(d);
        else que.push(pll(a, d));
    }
    vll ids(PA.size());
    iota(ALL(ids), 0);
    sort(ALL(ids), [&](int i, int j) { return -PD[i] > -PD[j]; });

    vll plus(M+1, -INF), minus(M+1, -INF);
    MonotoneConvexHullTrick<ll> cht;
    for (auto i : ids) cht.insert(-PD[i], PD[i] - PA[i] * 2);
    plus[0] = minus[0] = 0;
    for (int m = 1; m <= M; m++) {
        if (PA.size()) plus[m] = -cht.query_monotone(m) * m / 2;
        if (!que.empty()) {
            auto [s, d] = que.top();
            que.pop();
            minus[m] = minus[m-1] + s;
            que.push(pll(s+d, d));
        }
    }

    ll res = -INF;
    REP(m, M+1) chmax(res, plus[m] + minus[M-m]);
    cout << res << endl;
}
