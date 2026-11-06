#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// compression
template<class T> vector<T> compress(const vector<T> &v) {
    auto alt = v;
    sort(alt.begin(), alt.end());
    alt.erase(unique(alt.begin(), alt.end()), alt.end());
    vector<T> res(v.size());
    for (int i = 0; i < res.size(); i++) {
        res[i] = lower_bound(alt.begin(), alt.end(), v[i]) - alt.begin();
    }
    return res;
}
template<typename T, typename... Ts>
vector<T> merge_and_unique(const vector<T> &a, const Ts &...z) {
    vector<T> res = a;
    (res.insert(res.end(), z.begin(), z.end()), ...);
    sort(res.begin(), res.end());
    res.erase(unique(res.begin(), res.end()), res.end());
    return res;
}
