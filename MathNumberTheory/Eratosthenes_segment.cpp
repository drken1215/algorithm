//
// エラトステネスの区間篩
//
// cf.
//   高校数学の美しい物語: エラトスネテスの篩
//     https://mathtrain.jp/eratosthenes
//
// verified
//   AOJ 2858 Prime-Factor Prime
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2858
//


#include <bits/stdc++.h>
using namespace std;

const int MAX = 110000;
long long solve(long long L, long long R) {
    vector<bool> is_prime(MAX, true);
    is_prime[0] = is_prime[1] = false;
    vector<vector<long long>> prime_factors(R-L+1);
    for (long long p = 2; p * p <= R; ++p) {
        if (!is_prime[p]) continue;
        for (long long v = p*2; v < MAX; v += p) {
            is_prime[v] = false;
        }
        for (long long v = (L+p-1)/p*p; v <= R; v += p) {
            prime_factors[v-L].push_back(p);
        }
    }
    long long res = 0;
    for (long long v = L; v <= R; ++v) {
        long long prime_num = 0;
        long long v2 = v;
        const auto& ps = prime_factors[v2-L];
        for (auto p : ps) {
            while (v2 % p == 0) v2 /= p, ++prime_num;
        }
        if (v2 > 1) ++prime_num;
        if (is_prime[prime_num]) ++res;
    }
    return res;
}

int main() {
    long long L, R;
    cin >> L >> R;
    cout << solve(L, R) << endl;
}
