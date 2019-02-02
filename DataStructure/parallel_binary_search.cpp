//
// Partially Persistent Union-Find tree
//
// verified:
//   AGC 002 D - Stamp Rally
//     https://beta.atcoder.jp/contests/agc002/tasks/agc002_d
//


#include <iostream>
#include <vector>
using namespace std;

struct UnionFind {
    vector<int> par;
    
    UnionFind(int n) : par(n, -1) { }
    void init(int n) { par.assign(n, -1); }
    
    int root(int x) {
        if (par[x] < 0) return x;
        else return par[x] = root(par[x]);
    }
    
    bool issame(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y) {
        x = root(x); y = root(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y); // merge technique
        par[x] += par[y];
        par[y] = x;
        return true;
    }
    
    int size(int x) {
        return -par[root(x)];
    }
};

int main() {
    int N, M, Q;
    cin >> N >> M;
    vector<int> A(M), B(M);
    for (int i = 0; i < M; ++i) {
        scanf("%d %d", &A[i], &B[i]);
        --A[i], --B[i];
    }
    cin >> Q;
    vector<int> X(Q), Y(Q), Z(Q);
    for (int i = 0; i < Q; ++i) {
        scanf("%d %d %d", &X[i], &Y[i], &Z[i]);
        --X[i], --Y[i];
    }
    vector<int> le(Q), ri(Q);
    vector<vector<int> > vec(M);
    for (int i = 0; i < Q; ++i) le[i] = 0, ri[i] = M;
    while (true) {
        bool update = false;
        for (int i = 0; i < (int)vec.size(); ++i) vec[i].clear();
        for (int i = 0; i < Q; ++i) {
            if (ri[i] - le[i] > 1) {
                update = true;
                int mid = (le[i] + ri[i]) / 2;
                vec[mid].push_back(i);
            }
        }
        if (!update) break;
        UnionFind uf(N);
        for (int mid = 0; mid < M; ++mid) {
            for (auto q : vec[mid]) {
                int wa = uf.size(X[q]);
                if (!uf.issame(X[q], Y[q])) wa += uf.size(Y[q]);
                if (wa >= Z[q]) ri[q] = mid;
                else le[q] = mid;
            }
            uf.merge(A[mid], B[mid]);
        }
    }
    for (int q = 0; q < Q; ++q) cout << ri[q] << endl;
}
