#define PROBLEM "https://atcoder.jp/contests/abc480/tasks/abc480_g"

#include "sparse_table.hpp"
#include "monge_shortest_path_by_cht.hpp"

int main() {
    int N;
    cin >> N;
    long long INF = numeric_limits<long long>::max() / 2;
    long long M = numeric_limits<long long>::max() / (N * 3);
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    SparseTable<long long> stmi(A);
    SparseTable<long long> stma(A, [&](long long p, long long q){return max(p, q);});

    // Li Chao Tree は、通常、全定義域に anti-Monge 性を拡張する必要がある！！
    // もしくは、定義域を制限した区間挿入を使う
    auto func = [&](int l, int r) -> long long { 
        if (r < l + 2) return -M * ((l + 1 - r) * 2 + 1);
        else return stma.get(l, r) - stmi.get(l, r);
    };
    //assert(is_anti_monge(N + 1, N + 1, func));
    
    vector<int> start(N);
    for (int i = 0; i < N; i++) start[i] = i + 2;  // 辺 (i, i + 1) は禁止
    auto res = MongeShortestPathByCHT<long long>(N, func, start);
    cout << res.back().first << endl;
}
