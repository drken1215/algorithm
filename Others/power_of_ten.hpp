#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// power of ten
const vector<long long> TEN = {
    1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000,
    10000000000LL, 100000000000LL, 1000000000000LL, 10000000000000LL, 100000000000000LL,
    1000000000000000LL, 10000000000000000LL, 100000000000000000LL, 1000000000000000000LL
};
template<class T = long long> T get_ten(int n) {
    assert(n >= 0);
    if (n <= 18) return T(TEN[n]);
    T res = TEN[18];
    for (int i = 0; i < n - 18; i++) res *= 10;
    return res;
}
