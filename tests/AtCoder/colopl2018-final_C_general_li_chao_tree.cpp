#define PROBLEM "https://atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c"

#include <bits/stdc++.h>
using namespace std;

#include "general_li_chao_tree.hpp"

int main() {
    long long N;
    cin >> N;
    vector<long long> a(N), res(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    auto make_func = [&](long long i) {
        return [&, i](long long x) -> long long {
            return (x - i) * (x - i) + a[i];
        };
    };
    DynamicGeneralLiChaoTree<long long, decltype(make_func(0))> cht(0, 210000);
    for (long long i = 0; i < N; ++i) {
        cht.insert(make_func(i));
    }
    for (long long i = 0; i < N; ++i) res[i] = cht.query(i).first;
    for (int i = 0; i < N; ++i) cout << res[i] << endl;
}
