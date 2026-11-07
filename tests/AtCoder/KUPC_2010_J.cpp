#define PROBLEM "https://atcoder.jp/contests/kupc2012/tasks/kupc2012_10"

#include <bits/stdc++.h>
using namespace std;

#include "monge_knuth_yao.hpp"

int main() {
    int N;
    cin >> N;
    vector<long long> A(N), S(N+1, 0);
    for (int i = 0; i < N; i++) cin >> A[i], S[i+1] = S[i] + A[i];
    auto f = [&](int i, int j) -> long long { 
        if (j - i <= 1) return 0;
        else return S[j] - S[i];
    };
    auto res = KnuthYao<long long>(N, f);
    cout << res << endl;
}
