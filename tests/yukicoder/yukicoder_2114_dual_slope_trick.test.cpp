#define PROBLEM "https://yukicoder.me/problems/no/2114"

#include <bits/stdc++.h>
using namespace std;

#include "int128.hpp"
#include "dual_slope_trick.hpp"

// yukicoder No.2114 01 Matching
/*
    v1, ..., vN: 座標を小さい順に並べたもの。赤か青。赤を source、青を sink にする
    dp_{i}[x] := 最初の i 個を終えた時点で染み出しフローが x であるときの、染み出し分も含めたコストの最小値

    ・i 番目が赤のとき（D = x[i+1] - x[i] とする）
    nex[x] = dp[x-1] + D|x|（平行移動 + 直線加算）

    ・i 番目が青のとき（D = x[i+1] - x[i] とする）
    nex[x] = min(dp[x], dp[x+1]) + D|x|（スライド最小値 + 直線加算）
*/
int main() {
    const int RED = 0, BLUE = 1;
    long long N, M, K;
    cin >> N >> M >> K;
    vector<long long> B(N), R(M);
    map<long long, vector<pair<long long, int>>> mp;
    for (int i = 0; i < N; i++) cin >> B[i];
    for (int i = 0; i < M; i++) cin >> R[i];
    if (N < M) swap(N, M), swap(B, R);  // RED の方が少なくする
    for (int i = 0; i < N; i++) mp[B[i] % K].emplace_back(B[i] / K, BLUE);
    for (int i = 0; i < M; i++) mp[R[i] % K].emplace_back(R[i] / K, RED);
    
    auto calc = [&](vector<pair<long long, int>> &v) -> pair<bool, long long> {
        sort(v.begin(), v.end());
        DualSlopeTrick<long long, i128> st(1LL<<50);  // INF を積み上げても overflow しないように！
        for (int i = 0; i < v.size(); i++) {
            long long D = (i + 1 < v.size() ? v[i + 1].first - v[i].first : 0);
            if (v[i].second == RED) {
                st.slide(1);
                st.add_abs(D, 0);
            } else {
                st.slide(-1, 0);
                st.add_abs(D, 0);
            }
        }
        auto res = st.get_f0();
        return {bool(res < st.INF/2), res};
    };
    bool can = true;
    long long res = 0;
    for (auto [key, v] : mp) {
        auto [feasible, tmp] = calc(v);
        if (!feasible) can = false;
        else res += tmp;
    }
    cout << (can ? res : -1) << endl;
}
