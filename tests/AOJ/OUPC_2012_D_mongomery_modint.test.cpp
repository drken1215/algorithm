#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/problems/2353"

#include <bits/stdc++.h>
using namespace std;

#include "mongomery_modint.hpp"

int main() {
    using mint = MontgomeryModInt64;
    const long long MOD = 67280421310721LL;
    mint::set_mod(MOD);

    mint res = 0;
    int N;
    cin >> N;
    for (int i = 0; i < N; ++i) {
        long long o, y;
        cin >> o >> y;
        if (o == 1) res += y;
        else if (o == 2) res -= y;
        else if (o == 3) res *= y;
        else res /= y;
    }
    if (res.get() <= 1LL<<31) cout << res.get() << endl;
    else cout << (long long)res.get() - MOD << endl;
}
