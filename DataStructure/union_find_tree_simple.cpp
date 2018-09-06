//
// Union-Find tree
//
// verified:
//   AOJ Course DSL_1_A Set - Disjoint Set: Union Find Tree
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_1_A&lang=jp
//

// root(x): root node of x
// size(x): the size of the group of x


#include <iostream>
#include <vector>
#include <map>
using namespace std;

struct UnionFind {
    vector<int> par;
    
    UnionFind(int n) : par(n) { for (int i = 0; i < n; ++i) par[i] = i;  }
    void init(int n) { par.resize(n); for (int i = 0; i < n; ++i) par[i] = i; }
    
    int root(int x) {
        if (par[x] == x) return x;
        else return par[x] = root(par[x]);
    }
    
    bool issame(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y) {
        x = root(x); y = root(y);
        if (x == y) return false;
        par[x] = y;
        return true;
    }
};



int main() {
    int N, Q; cin >> N >> Q;
    UnionFind uf(N);
    for (int q = 0; q < Q; ++q) {
        int com, x, y; cin >> com >> x >> y;
        if (com == 0) uf.merge(x, y);
        else {
            if (uf.issame(x, y)) cout << 1 << endl;
            else cout << 0 << endl;
        }
    }
}
