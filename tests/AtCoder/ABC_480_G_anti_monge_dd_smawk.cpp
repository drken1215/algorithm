#define PROBLEM "https://atcoder.jp/contests/abc480/tasks/abc480_g"

#include "sparse_table.hpp"
#include "anti_monge_dd_smawk.hpp"

int main() {
    int N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    SparseTable<long long> stmi(A);
    SparseTable<long long> stma(A, [&](long long p, long long q){return max(p, q);});

    auto func = [&](int l, int r) -> long long {
        assert(r - l >= 2);
        return stma.get(l, r) - stmi.get(l, r);
    };
    auto res = antiMongeDDSMAWK<long long>(N, func, 2);
    cout << res.back().first << endl;
}
