#define PROBLEM "https://judge.yosupo.jp/problem/many_aplusb_128bit"

#include <bits/stdc++.h>
using namespace std;

#include "int128.hpp"
#include "fast_io.hpp"

int main() {
    FastRead Read;
    FastWrite Write;
    
    int T;
    Read(T);
    for (int t = 0; t < T; ++t) {
        i128 a, b;
        Read(a, b);
        Write.ln(a + b);
    }
}
