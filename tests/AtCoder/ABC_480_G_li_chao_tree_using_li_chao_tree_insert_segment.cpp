#define PROBLEM "https://atcoder.jp/contests/abc480/tasks/abc480_g"

#include "template.hpp"
#include "sparse_table.hpp"
#include "general_li_chao_tree.hpp"

// Li-Chao Tree の「区間挿入」を用いる O(N log^2 N) 解法
int main() {
    int N;
    cin >> N;
    vll A(N);
    cin >> A;

    SparseTable<ll> st1(A, [&](ll p, ll q){return max(p, q);}), st2(A);
    auto func = [&](int l, int r) -> ll {
        return st1.get(l, r) - st2.get(l, r);
    };

    // dp[j] = min_{j < i-1} (dp[i] + func(i, j))
    // f_{i}(x) := dp[i] + func(i, x) として。定義域を x >= i + 2 とする
    vector<long long> dp(N + 1);
    auto make_func = [&](int i) {
        return [&, i](int x) -> ll {
            return dp[i] + func(i, x);
        };
    };
    DynamicGeneralLiChaoTree<ll, decltype(make_func(0))> cht(0, N + 1);
    dp[0] = 0;
    cht.insert_segment(make_func(0), 2, N + 1);
    for (int i = 2; i <= N; i++) {
        dp[i] = cht.query(i).first;
        if (i + 2 <= N) cht.insert_segment(make_func(i), i + 2, N + 1);
    }
    cout << dp[N] << endl;
}
