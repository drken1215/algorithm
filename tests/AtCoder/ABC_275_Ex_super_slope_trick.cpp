#define PROBLEM "https://atcoder.jp/contests/abc275/tasks/abc275_h"

#include <bits/stdc++.h>
using namespace std;

#include "super_slope_trick.hpp"

// ABC 275 Ex - Monster
/*
　木を Cartesian 木（根 = 区間の最大 B）で辿り、各ノードで
　dp[v][x] := v を根とする部分木に属するモンスターについて、
             残り x 以下の攻撃を要する状態にするための消費魔力の総和の最小値

　更新時
　dp[l], dp[r] のうち x < A[v] の領域を削る
　dp[v][x] = min_{y >= x} (dp[l][y] + dp[r][y] + B[v]y) - B[v]x
　その後、x ≧ A[v] に制限された定義域を、x ≧ 0 まで広げる

　答え
　eval(dp[root], 0)
*/
int main() {
    long long N;
    cin >> N;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> B[i];
    long long MAXX = 0;
    for (auto a : A) MAXX = max(MAXX, a);

    // Cartesian 木の構築とトポロジカルソート
    vector<int> lch(N, -1), rch(N, -1), stk;
    for (int i = 0; i < N; i++) {
        int last = -1;
        while (!stk.empty() && B[stk.back()] < B[i]) { last = stk.back(); stk.pop_back(); }
        lch[i] = last;
        if (!stk.empty()) rch[stk.back()] = i;
        stk.push_back(i);
    }
    int root = stk.front();
    vector<int> order;
    vector<int> s1 = {root};
    while (!s1.empty()) {
        int v = s1.back();
        s1.pop_back();
        order.push_back(v);
        if (lch[v] != -1) s1.push_back(lch[v]);
        if (rch[v] != -1) s1.push_back(rch[v]);
    }
    reverse(order.begin(), order.end());

    Slope_Trick_Super<long long> ST;
    using F = typename decltype(ST)::FUNC;
    vector<F> dp(N);
    for (int v : order) {
        // 定義域 [0, MAXX] 上の f(x) 0
        dp[v] = ST.segment_func(0, MAXX, 0, 0);

        // 左右からのマージ
        if (lch[v] != -1) dp[v] = ST.add(dp[v], dp[lch[v]]);
        if (rch[v] != -1) dp[v] = ST.add(dp[v], dp[rch[v]]);

        // 定義域を [A[v], MAXX] に制限
        dp[v] = ST.restrict_domain(dp[v], A[v], MAXX);

        // その後の処理
        dp[v] = ST.add_linear(dp[v], B[v], 0);
        dp[v] = ST.clear_left(dp[v]);
        dp[v] = ST.add_linear(dp[v], -B[v], 0);

        // 定義域を傾き -B[v] で [0, MAXX] に拡張
        dp[v] = ST.extend_left(dp[v], 0, -B[v]);

        //COUT(v); ST.dump(dp[v]);
    }
    cout << ST.eval(dp[root], 0) << "\n";                       // 答えは eval(・, 0)
}
