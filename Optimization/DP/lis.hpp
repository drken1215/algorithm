#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// dp[i] := 長さが i の増加部分列として最後尾の要素のとりうる最小値
template<class T> int find_lis(vector<T> a,  bool is_strong = true) {
    const T INF = numeric_limits<T>::max() / 2;
    int n = (int)a.size();
    vector<T> dp(n, INF);
    for (int i = 0; i < n; ++i) {
        if (is_strong) *lower_bound(dp.begin(), dp.end(), a[i]) = a[i];
        else *upper_bound(dp.begin(), dp.end(), a[i]) = a[i];
    }
    return lower_bound(dp.begin(), dp.end(), INF) - dp.begin();
}
