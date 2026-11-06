#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// next combination
template<class T> bool next_combination(T &bit, int N) {
    T x = bit & -bit, y = bit + x;
    bit = (((bit & ~y) / x) >> 1) | y;
    return (bit < (1LL << N));
}
