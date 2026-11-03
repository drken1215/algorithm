#define PROBLEM "https://atcoder.jp/contests/abc334/tasks/abc334_b"

#include <bits/stdc++.h>
using namespace std;

#include "floor_ceil.hpp"

int main() {
    long long A, M, L, R;
    cin >> A >> M >> L >> R;
    long long l = ceil(L - A, M), r = floor(R - A, M);
    cout << r - l + 1 << endl;
}