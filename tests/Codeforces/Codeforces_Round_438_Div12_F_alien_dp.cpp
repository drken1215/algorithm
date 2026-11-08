#define PROBLEM "https://codeforces.com/contest/868/problem/F"

#include <bits/stdc++.h>
using namespace std;

#include "alien_dp.hpp"

int main() {
    long long N, K;
    cin >> N >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];

    struct STATE {
        vector<long long> nums;
        long long cost;
        STATE(int N) : nums(N + 1, 0), cost(0) {}
    };
    STATE ini(N);
    auto add = [&](STATE &st, int i) -> void {
        st.cost += st.nums[A[i]];
        st.nums[A[i]]++;
    };
    auto del = [&](STATE &st, int i) -> void {
        st.nums[A[i]]--;
        st.cost -= st.nums[A[i]];
    };
    auto get = [&](STATE &st) -> long long {
        return st.cost;
    };

    AlienDP<long long> adp;
    auto res = adp.solve(N, ini, add, del, get, K, "simple_larsch");  // we should use slide access
    cout << res.first << endl;
}
