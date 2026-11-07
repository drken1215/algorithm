#define PROBLEM "https://atcoder.jp/contests/abc224/tasks/abc224_b"

#include <bits/stdc++.h>
using namespace std;

#include "stl_container_stream.hpp"
#include "monge_checker.hpp"

int main() {
    int H, W;
    cin >> H >> W;
    vector A(H, vector(W, 0LL));
    cin >> A;
    auto f = [&](int i, int j) -> long long { return A[i][j]; };
    Yes(is_monge(H, W, f, false, false));
}
