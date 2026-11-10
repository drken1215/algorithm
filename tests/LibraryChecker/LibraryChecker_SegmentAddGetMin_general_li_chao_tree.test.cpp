#define PROBLEM "https://judge.yosupo.jp/problem/segment_add_get_min"

#include <bits/stdc++.h>
using namespace std;

#include "fast_io.hpp"
#include "general_li_chao_tree.hpp"

int main() {
    FastRead Read;
    FastWrite Write;

    // クエリ先読み (座標圧縮のため)
    long long N, Q;
    Read(N, Q);
    vector<long long> typ(Q), l(N + Q), r(N + Q), a(N + Q), b(N + Q), p(Q), xs;
    for (int i = 0; i < N; i++) {
        Read(l[i], r[i], a[i], b[i]);
    }
    for (int q = 0; q < Q; q++) {
        Read(typ[q]);
        if (typ[q] == 0) {
            Read(l[q + N], r[q + N], a[q + N], b[q + N]);
        } else {
            Read(p[q]);
            xs.emplace_back(p[q]);
        }
    }
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());

    // 関数群を定義
    auto make_func = [&](long long a, long long b) {
        return [&, a, b](long long x) -> long long {
            return a * x + b;
        };
    };

    // 座標圧縮版 Li-Chao Tree
    GeneralLiChaoTree<long long, decltype(make_func(0, 0))> cht(xs);
    for (int i = 0; i < N; i++) {
        cht.insert_segment(make_func(a[i], b[i]), l[i], r[i]);
    }
    for (int q = 0; q < Q; q++) {
        if (typ[q] == 0) {
            cht.insert_segment(make_func(a[q + N], b[q + N]), l[q + N], r[q + N]);
        } else {
            auto [res, id] = cht.query(p[q]);
            if (id == -1) Write.ln("INFINITY");
            else Write.ln(res);
        }
    }
}
