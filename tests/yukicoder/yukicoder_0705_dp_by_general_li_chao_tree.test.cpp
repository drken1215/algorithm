#define PROBLEM "https://yukicoder.me/problems/no/705"

#include <bits/stdc++.h>
using namespace std;

#include "general_li_chao_tree.hpp"

// yukicoder No.705 ゴミ拾い Hard
/*
    A, X, Y: N 個
    これらを区間に分割していく
    　dp[j] = min_{0 ≦ i < j}(dp[i] + |A[j-1] - X[i]|^3 + |-Y[i]|^3)
    i -> j のコスト：|A[j-1] - X[i]|^3 + |-Y[i]|^3 ...... 差の凸関数 (Monge) + 縞々 (Monge) -> Monge
    スタート: 0, ゴール: N
*/
int main() {
    int N;
    cin >> N;
    vector<long long> A(N), X(N), Y(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> X[i];
    for (int i = 0; i < N; i++) cin >> Y[i];

    vector<long long> dp(N);
    auto func = [&](int i, int j) -> long long {
        long long dx = abs(A[j-1] - X[i]), dy = abs(Y[i]);
        return dx * dx * dx + dy * dy * dy;
    };
    auto make_func = [&](int i) {
        return [&, i](int j) -> long long {
            return dp[i] + func(i, j);
        };
    };
    GeneralLiChaoTree<long long, decltype(make_func(0))> cht(0, N + 1);
    dp[0] = 0;
    for (int j = 1; j <= N; j++) {
        cht.insert_segment(make_func(j - 1), j, N + 1);
        dp[j] = cht.query(j).first;
    }
    cout << dp[N] << endl;
}
