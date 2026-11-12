#define PROBLEM "https://yukicoder.me/problems/no/2546"

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "general_li_chao_tree.hpp"

/*
    f[a,d](m) = (A*2 + D(m-1)) * m / 2
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

    vll plus(M+1, -INF), minus(M+1, -INF);
    auto make_func = [&](ll a, ll d) {
        return [&, a, d](ll m) -> ll {
            return (-d * (m - 1) - a * 2) * m / 2;
        };
    };
    GeneralLiChaoTree<ll, decltype(make_func(0, 0))> cht(0, M+1);
    REP(i, PA.size()) cht.insert(make_func(PA[i], PD[i]));
    plus[0] = minus[0] = 0;
    for (int m = 1; m <= M; m++) {
        if (PA.size()) plus[m] = -cht.query(m).first;
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
