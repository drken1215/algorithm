#define PROBLEM "https://yukicoder.me/problems/no/2404"

#include <bits/stdc++.h>
using namespace std;

#include "template.hpp"
#include "general_li_chao_tree.hpp"

/*
    h(x) = a(x - s)(x - t)
*/
int main() {
    long long a, Q;
    cin >> a >> Q;
    auto make_func = [&](ll p, ll q) {
        return [&, p, q](ll t) -> ll {
            return a * (t - p) * (t - q);
        };
    };
    DynamicGeneralLiChaoTree<ll, decltype(make_func(0, 0))> cht(0, 123456789);
    REP(q, Q) {
        ll typ, s, t; cin >> typ;
        if (typ == 1) {
            cin >> s >> t;
            cht.insert(make_func(s, t));
        } else {
            cin >> t;
            cout << max(-cht.query(t).first, 0LL) << endl;
        }
    }
}
