#define PROBLEM "https://atcoder.jp/contests/abc407/tasks/abc407_f"

#include "sparse_table.hpp"

int main() {
    int N;
    cin >> N;
    using Node = pair<long long, int>;
    vector<Node> A(N);
    for (int i = 0; i < N; i++) {
        cin >> A[i].first;
        A[i].second = i;
    }
    SparseTable<Node> st(A, [&](Node a, Node b){return max(a, b);});

    // imos method
    vector<long long> res(N+3, 0);
    auto imos = [&](int l, int r, long long fval, long long d) -> void {
        // 区間 [l, r) に、初項 fval, 公差 d の等差数列を足す
        long long lval = fval + d * (r - l - 1);
        res[l] += fval, res[l+1] += -fval + d, res[r] += -lval - d, res[r+1] += lval;
    };

    // 再帰的に最大値を足し込む
    auto rec = [&](auto &&rec, int l, int r) -> void {
        if (r - l <= 0) return;
        auto [val, m] = st.get(l, r);
        int L = min(m+1-l, r-m), R = max(m+1-l, r-m);
        imos(1, L+1, val, val);
        imos(L+1, R+1, val * L, 0);
        imos(R+1, r-l+1, val * (L-1), -val);
        rec(rec, l, m);
        rec(rec, m+1, r);
    };
    rec(rec, 0, N);

    // imos の復元 (2 回累積和をとる)
    for (int i = 0; i+1 < res.size(); i++) res[i+1] += res[i];
    for (int i = 0; i+1 < res.size(); i++) res[i+1] += res[i];

    // 出力
    for (int k = 1; k <= N; k++) cout << res[k] << '\n';
}
