#define PROBLEM "https://atcoder.jp/contests/typical90/tasks/typical90_a"

#include <bits/stdc++.h>
using namespace std;

#include "binary_search.hpp"

int main() {
    long long N, L, K;
    cin >> N >> L >> K;
    vector<long long> A(N+2, 0);
    for (int i = 0; i < N; i++) cin >> A[i+1];
    A[N+1] = L;

    // 長さ x 以上ずつ切って、K+1 個を切り出せるか
    auto check = [&](long long x) -> bool {
        int prev = 0;
        long long num = 0;
        for (int i = 0; i < N+2; i++) {
            if (A[i] - A[prev] >= x) {
                prev = i;
                num++;
            }
        }
        return num >= K+1;
    };
    long long res = binary_search_10(1LL, L, check);
    cout << res << endl;
}
