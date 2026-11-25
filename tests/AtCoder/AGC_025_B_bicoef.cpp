#define PROBLEM "https://atcoder.jp/contests/agc025/tasks/agc025_b"

#include <bits/stdc++.h>
using namespace std;

#include "bicoef.hpp"

int main() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;

    long long N, A, B, K;
    cin >> N >> A >> B >> K;

    BiCoef<mint> bc(N + 1);
    mint res = 0;
    for (long long a = 0; a <= N; ++a) { // A の個数を全探索
        long long rem = K - a * A;
        if (rem % B != 0) continue;
        long long b = rem / B;
        if (b > N) continue;
        mint tmp = bc.com(N, a) * bc.com(N, b);
        res += tmp;
    }   
    cout << res << endl; 
}
