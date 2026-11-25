#define PROBLEM "https://judge.yosupo.jp/problem/binomial_coefficient_prime_mod"

#include <bits/stdc++.h>
using namespace std;

#include "dynamic_modint.hpp"
#include "bicoef.hpp"

int main() {
    const int MAX = 11000000;
    int T, M;
    cin >> T >> M;
    using mint = DynamicModint;
    mint::set_mod(M);
    
    BiCoef<mint> bc(MAX);
    while (T--) {
        int n, k;
        cin >> n >> k;
        cout << bc.com(n, k) << endl;
    }
}
