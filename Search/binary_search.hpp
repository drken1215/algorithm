#pragma once
#include <bits/stdc++.h>
using namespace std;

// binary search
template<class T, class FUNC> T binary_search_01(T ng, T ok, const FUNC &f) {
    while (abs(ok - ng) > 1) {
        T mid = ok + ng >> 1;
        if (f(mid)) ok = mid;
        else ng = mid;
    }
    return ok;
}
template<class T, class FUNC> T binary_search_10(T ok, T ng, const FUNC &f) {
    return binary_search_01(ng, ok, f);
}
