//
// Monge 単一始点最短路 by LARSCH 法 in O(N)
//   頂点数 N+1 の DAG, 頂点 i, j 間のコスト f(i, j) が Monge であることを仮定
//
// verified
//   AtCoder EDPC Z - Frog 3
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//
//   Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
//     https://codeforces.com/contest/319/problem/C 
// 
//   yukicoder No.705 ゴミ拾い Hard
//     https://yukicoder.me/problems/no/705 
//
// Reference:
//   noshi: 簡易版 LARSCH Algorithm
//     https://noshi91.hatenablog.com/entry/2023/02/18/005856
//  


#include <bits/stdc++.h>
using namespace std;


// LARSCH, in O(N)
// vertex: 0, 1, 2, ..., N, f(i, j) must be Monge
template<class VAL> struct LARSCH {

};

template<class VAL> struct MongeShortestPath {
    VAL INF = numeric_limits<VAL>::max() / 2;
    int CNT_INF = numeric_limits<int>::max() / 2;

    // results
    vector<VAL> dp;
    vector<int> cnt, prev;

    // solver (random access ver)
    template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f, bool minimize_cnt = true) {
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;

        auto relax = [&](int l, int r) -> void {
            VAL val = dp[l] + f(l, r);
            int cn = cnt[l] + 1;
            if (dp[r] > val || (dp[r] == val && (minimize_cnt ? cn < cnt[r] : cn > cnt[r]))) {
                dp[r] = val;
                cnt[r] = cn;
                prev[r] = l;
            }
        };
        auto rec = [&](auto &&rec, int l, int r) -> void {
            if (r - l <= 1) return;
            int m = (l + r) / 2;
            for (int k = prev[l]; k <= prev[r]; k++) relax(k, m);
            rec(rec, l, m);
            for (int k = l + 1; k <= m; k++) relax(k, r);
            rec(rec, m, r);
        };

        if (N > 0) relax(0, N), rec(rec, 0, N);
        vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
        res[0].first = VAL(0);
        for (int i = 1; i <= N; i++) res[i] = {dp[i], prev[i]};
        return res;
    }

    // solver (slide access ver)
    template<class STATE, class ADD, class DEL, class GETCOST> vector<pair<VAL, int>> solve
    (int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, bool minimize_cnt = true) {
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;

        STATE win[2] = { ini, ini };
        int cl[2] = {0, 0}, cr[2] = {0, 0};

        auto move_cursor = [&](int c, int l, int r) {
            while (cr[c] < r) add(win[c], cr[c]++);
            while (cl[c] > l) add(win[c], --cl[c]);
            while (cr[c] > r) del(win[c], --cr[c]);
            while (cl[c] < l) del(win[c], cl[c]++);
        };
        auto relax = [&](int c, int l, int r) {
            move_cursor(c, l, r);
            VAL val = dp[l] + get(win[c]);
            int cn = cnt[l] + 1;
            if (dp[r] > val || (dp[r] == val && (minimize_cnt ? cn < cnt[r] : cn > cnt[r]))) {
                dp[r] = val;
                cnt[r] = cn;
                prev[r] = l;
            }
        };
        auto rec = [&](auto &&rec, int l, int r) -> void {
            if (r - l <= 1) return;
            int m = (l + r) / 2;
            for (int k = prev[l]; k <= prev[r]; k++) relax(0, k, m);
            rec(rec, l, m);
            for (int k = l + 1; k <= m; k++) relax(1, k, r);
            rec(rec, m, r);
        };

        if (N > 0) relax(0, 0, N), rec(rec, 0, N);
        vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
        res[0].first = VAL(0);
        for (int i = 1; i <= N; i++) res[i] = {dp[i], prev[i]};
        return res;
    }

    // 遷移回数を問わない場合
// template <typename T, typename F>
// vc<T> monge_shortest_path(int N, F f) {
//   vc<T> dp(N + 1, infty<T>);
//   dp[0] = 0;
//   LARSCH<T> larsch(N, [&](int i, int j) -> T {
//     ++i;
//     if (i <= j) return infty<T>;
//     return dp[j] + f(j, i);
//   });
//   FOR(r, 1, N + 1) {
//     int l = larsch.get_argmin();
//     dp[r] = dp[l] + f(l, r);
//   }
//   return dp;
// }

    vector<int> reconstruct() {
        int N = (int)dp.size() - 1;
        vector<int> path;
        for (int v = N; v > 0; v = prev[v]) path.emplace_back(v);
        path.emplace_back(0);
        reverse(path.begin(), path.end());
        return path;
    }
};


//------------------------------//
// Examples
//------------------------------//

// AtCoder EDPC Z - Frog 3
/*
    H は単調増加数列
    chmin(dp[j], dp[i] + (H[j] - H[i])^2 + C)
    i -> j のコスト：(H[j] - H[i])^2 ...... 差の凸関数は Monge
    スタート: 0, ゴール: N-1
*/
void EDPC_Z() {
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N);
    for (long long i = 0; i < N; i++) cin >> H[i];
    auto func = [&](int i, int j) -> long long {
        return (H[j] - H[i]) * (H[j] - H[i]) + C;
    };
    MongeShortestPath<long long> msp;
    auto res = msp.solve(N-1, func);
    cout << res[N-1].first << endl;
}

// Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
/*
    A: 単調増加, B: 単調減少, ともに長さ N
    i -> j のコストが、B[i] × A[j] で与えられる ..... 単調増加 × 単調減少は Monge
    スタート: 0, ゴール: N-1
*/
void Codeforces_189_C() {
    long long N;
    cin >> N;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> B[i];
    auto func = [&](int i, int j) -> long long {
        return B[i] * A[j];
    };
    MongeShortestPath<long long> msp;
    auto res = msp.solve(N-1, func);
    cout << res[N-1].first << endl;
}

// yukicoder No.705 ゴミ拾い Hard
/*
    A, X, Y: N 個
    これらを区間に分割していく
    　dp[j] = min_{0 ≦ i < j}(dp[i] + |A[j-1] - X[i]|^3 + |-Y[i]|^3)
    i -> j のコスト：|A[j-1] - X[i]|^3 + |-Y[i]|^3 ...... 差の凸関数 (Monge) + 縞々 (Monge) -> Monge
    スタート: 0, ゴール: N
*/
void yukicoder_705() {
    int N;
    cin >> N;
    vector<long long> A(N), X(N), Y(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> X[i];
    for (int i = 0; i < N; i++) cin >> Y[i];
    auto func = [&](int i, int j) -> long long {
        long long dx = abs(A[j-1] - X[i]), dy = abs(Y[i]);
        return dx * dx * dx + dy * dy * dy;
    };
    MongeShortestPath<long long> msp;
    auto res = msp.solve(N, func);
    cout << res[N].first << endl;
}


int main() {
    EDPC_Z();
    //Codeforces_189_C();
    //yukicoder_705();
}