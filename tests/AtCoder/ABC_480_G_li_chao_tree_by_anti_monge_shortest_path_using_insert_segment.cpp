#define PROBLEM "https://atcoder.jp/contests/abc480/tasks/abc480_g"

#include "sparse_table.hpp"
#include "monge_shortest_path_by_cht.hpp"

// 内部で Li-Chao Tree の「区間挿入」を用いる O(N log^2 N) 解法
int main() {
    int N;
    cin >> N;
    long long INF = numeric_limits<long long>::max() / 2;
    long long M = numeric_limits<long long>::max() / (N * 3);
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    SparseTable<long long> stmi(A);
    SparseTable<long long> stma(A, [&](long long p, long long q){return max(p, q);});
    
    vector<int> start(N);
    for (int i = 0; i < N; i++) start[i] = i + 2;  // 辺 (i, i + 1) は禁止

    // Li Chao Tree は、通常、全定義域に anti-Monge 性を拡張する必要がある！！
    // もしくは、定義域を制限した区間挿入を使う
    auto func = [&](int l, int r) -> long long { 
        assert(r >= start[l]);
        return stma.get(l, r) - stmi.get(l, r);
    };
    auto res = MongeShortestPathByCHT<long long>(N, func, start, true);
    cout << res.back().first << endl;
}
