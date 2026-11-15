#define PROBLEM "https://atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c"

#include <bits/stdc++.h>
using namespace std;

#include "monotone_convex_hull_trick.hpp"

int main() {
    long long N;
    cin >> N;
    vector<long long> a(N), res(N, 1LL<<60);
    for (int i = 0; i < N; ++i) cin >> a[i];
    MonotoneConvexHullTrick<long long> cht;
    for (long long i = 0; i < N; ++i) cht.insert(-2LL * i, a[i] + i*i);
    for (long long i = 0; i < N; ++i) res[i] = cht.query_monotone(i) + i*i;
    for (int i = 0; i < N; ++i) cout << res[i] << endl;
}
