#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// prev subset
template<class T> bool prev_subset(T &S, T U) {
    assert(S >= 0 && U >= 0);
    if (S == 0) return false;
    S = (S - 1) & U;
    return true;
}
