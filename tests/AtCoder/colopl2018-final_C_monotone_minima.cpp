#define PROBLEM "https://atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c"

#include <bits/stdc++.h>
using namespace std;

#include "monotone_minima.hpp"

// COLOCON 2018 Final C - スペースエクスプローラー高橋君
/*
    M[i][j] = (j - i)^2 は Monge
    列 j に一様に A[j] を加えても Monge
*/
int main() {
    long long N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    auto func = [&](long long i, long long j) -> long long {
        return A[j] + (j - i) * (j - i);
    };
    auto res = MonotoneMinima<long long>(N, N, func);
    for (int i = 0; i < N; i++) cout << res[i].first << endl;
}
