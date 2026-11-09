#define PROBLEM "https://judge.yosupo.jp/problem/segment_add_get_min"

#include <bits/stdc++.h>
using namespace std;

#include "fast_io.hpp"
#include "general_li_chao_tree.hpp"

int main() {
    FastRead Read;
    FastWrite Write;

    long long LX = -(1LL << 30), RX = 1LL << 30;
    auto make_func = [&](long long a, long long b) {
        return [&, a, b](long long x) -> long long {
            return a * x + b;
        };
    };
    DynamicGeneralLiChaoTree<long long, decltype(make_func(0, 0))> cht(LX, RX);
    long long N, Q, typ, l, r, a, b, p;
    Read(N, Q);
    for (int i = 0; i < N; i++) {
        Read(l, r, a, b);
        cht.insert_segment(make_func(a, b), l, r);
    }
    while (Q--) {
        Read(typ);
        if (typ == 0) {
            Read(l, r, a, b);
            cht.insert_segment(make_func(a, b), l, r);
        } else {
            Read(p);
            auto [res, id] = cht.query(p);
            if (id == -1) Write.ln("INFINITY");
            else Write.ln(res);
        }
    }
}
