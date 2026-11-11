#define PROBLEM "https://yukicoder.me/problems/no/409"

#include <bits/stdc++.h>
using namespace std;

#include "macro.hpp"
#include "stl_container_stream.hpp"
#include "monge_shortest_path_by_cht.hpp"

int main() {
    ll N, A, B, W;
    cin >> N >> A >> B >> W;
    vll D(N); 
    cin >> D;
    auto func = [&](ll i, ll j) -> ll {
        ll len = j - i;
        return D[j-1] - A * (len - 1) + B * len * (len - 1) / 2;
    };
    auto x = MongeShortestPathByCHT<long long>(N, func);
    ll res = x.back().first + W;
    REP(i, N+1) {
        ll w = x[i].first - A * (N - i) + B * (N - i) * (N - i + 1) / 2 + W;
        chmin(res, w);
    }
    cout << res << endl;
}