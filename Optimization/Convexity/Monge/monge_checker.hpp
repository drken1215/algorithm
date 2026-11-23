#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// check whether N x M matrix f is Monge or not, in O(NM)
template<class FUNC> bool is_monge
(int N, int M, const FUNC &f, bool upper_triangular = false, bool output_detail = true) {
    assert(N > 0 && M > 0);
    for (int i = 0; i + 1 < N; i++) {
        for (int j = (upper_triangular ? i+2 : 0); j + 1 < M; j++) {
            auto f00 = f(i, j), f01 = f(i, j + 1), f10 = f(i + 1, j), f11 = f(i + 1, j + 1);
            if (f00 + f11 > f10 + f01) {
                if (output_detail) {
                    cout << "Not Monge! " << endl;
                    cout << "f(" << i << ", " << j << ") = " << f(i, j) << ", ";
                    cout << "f(" << i << ", " << j + 1 << ") = " << f(i, j + 1) << endl;
                    cout << "f(" << i + 1 << ", " << j << ") = " << f(i + 1, j) << ", ";
                    cout << "f(" << i + 1 << ", " << j + 1 << ") = " << f(i + 1, j + 1) << endl;
                }
                return false;
            }
        }
    }
    return true;
}

// check whether N x M matrix f is anti-Monge or not, in O(NM)
template<class FUNC> bool is_anti_monge
(int N, int M, const FUNC &f, bool upper_triangular = false, bool output_detail = true) {
    assert(N > 0 && M > 0);
    for (int i = 0; i + 1 < N; i++) {
        for (int j = (upper_triangular ? i+2 : 0); j + 1 < M; j++) {
            auto f00 = f(i, j), f01 = f(i, j + 1), f10 = f(i + 1, j), f11 = f(i + 1, j + 1);
            if (f00 + f11 < f10 + f01) {
                if (output_detail) {
                    cout << "Not anti-Monge! " << endl;
                    cout << "f(" << i << ", " << j << ") = " << f(i, j) << ", ";
                    cout << "f(" << i << ", " << j + 1 << ") = " << f(i, j + 1) << endl;
                    cout << "f(" << i + 1 << ", " << j << ") = " << f(i + 1, j) << ", ";
                    cout << "f(" << i + 1 << ", " << j + 1 << ") = " << f(i + 1, j + 1) << endl;
                }
                return false;
            }
        }
    }
    return true;
}
