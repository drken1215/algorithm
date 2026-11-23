#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// xorshift による乱数生成、周期は2^128-1
struct xorshift {
    unsigned int tx, ty, tz, tw;
    xorshift() : tx(123456789), ty(362436069), tz(521288629), tw(88675123) {}

    // integer
    unsigned int rand_int() {
        unsigned int tt = tx ^ (tx << 11);
        tx = ty, ty = tz, tz = tw;
        return tw = (tw ^ (tw >> 19)) ^ (tt ^ (tt >> 8));
    }
    int rand_int(int left, int right) {
        assert(left < right);
        return rand_int() % (right - left) + left;
    }
    long long rand_ll(long long left, long long right) {
        long long a = rand_int(), b = rand_int();
        return (a * numeric_limits<unsigned int>::max() + b) % (right - left) + left;
    }

    // double
    double rand_double() {
        return (double)rand_int() / numeric_limits<unsigned int>::max();
    }

    // random shffle
    template<class T> void shuffle(vector<T> &vec) {
        int n = (int)vec.size();
        for (int i = n - 1; i > 0; i--) {
            int k = rand_int() % (i + 1);
            swap(vec[i], vec[k]);
        }
    }
};
