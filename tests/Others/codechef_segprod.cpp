#define PROBLEM "https://www.codechef.com/problems/SEGPROD"

#include "disjoint_sparse_table.hpp"

int main() {
    int T;
    scanf("%d", &T);
    int N, P, Q;
    DisjointSparseTable<long long> dst([&](long long a, long long b){return a*b%P;});
    for (int CASE = 0; CASE < T; ++CASE) {
        scanf("%d %d %d", &N, &P, &Q);
        vector<long long> v(N);
        for (int i = 0; i < N; ++i) scanf("%lld", &v[i]);
        dst.init(v);
        vector<int> b(Q/64+2);
        for (int i = 0; i < (int)b.size(); ++i) scanf("%d", &b[i]);
        int x = 0, L = 0, R = 0;
        for (int query = 0; query < Q; ++query) {
            if (query % 64 == 0) L = (b[query/64] + x) % N, R = (b[query/64+1] + x) % N;
            else L = (L + x) % N, R = (R + x) % N;
            if (L > R) swap(L, R);
            x = (dst.get(L, R+1) + 1) % P;
        }
        printf("%d\n", x);
    }
}
