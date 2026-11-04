#define PROBLEM "https://judge.yosupo.jp/problem/line_add_get_min"

#include <bits/stdc++.h>
using namespace std;

#include "general_li_chao_tree.hpp"

int main() {
    long long LX = -(1LL << 30), RX = 1LL << 30;
    auto make_func = [&](long long a, long long b) {
        return [&, a, b](long long x) -> long long {
            return a * x + b;
        };
    };
    DynamicGeneralLiChaoTree<long long, decltype(make_func(0, 0))> cht(LX, RX);
    long long N, Q, typ, a, b, p;
    cin >> N >> Q;
    for (int i = 0; i < N; i++) {
        cin >> a >> b;
        cht.insert(make_func(a, b));
    }
    while (Q--) {
        cin >> typ;
        if (typ == 0) {
            cin >> a >> b;
            cht.insert(make_func(a, b));
        } else {
            cin >> p;
            auto [res, id] = cht.query(p);
            cout << res << endl;
        }
    }
}
