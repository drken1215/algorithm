#define PROBLEM "https://atcoder.jp/contests/abc127/tasks/abc127_e"

#include <bits/stdc++.h>
using namespace std;

#include "bicoef.hpp"

int main() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;

    long long N, M, K;
    cin >> N >> M >> K;

    BiCoef<mint> bc(N * M);
    mint sum = 0;
    for (int i = 0; i <= N-1; ++i) {
        for (int j = 0; j <= M-1; ++j) {
            mint tmp = mint(N - i) * mint(M - j) * mint(i + j);
            if (i != 0 && j != 0) tmp *= 2;
            sum += tmp;
        }
    }
    cout << sum * bc.com(N * M - 2, K - 2) << endl;
}
