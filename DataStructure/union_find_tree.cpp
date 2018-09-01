//
// Projction and Reflection
//
// verified:
//   AtCoder Typical Contest 001 B - Union Find
//

// root(x): root node of x
// size(x): the size of the group of x

#include <iostream>
#include <vector>
#include <map>
using namespace std;

struct UnionFind {
    vector<int> par, rank, sz;
    
    UnionFind(int n) : par(n), rank(n, 0), sz(n, 1) {
        for (int i = 0; i < n; ++i) par[i] = i;
    }
    void init(int n) {
        par.resize(n); rank.resize(n); sz.resize(n);
        for (int i = 0; i < n; ++i) par[i] = i, rank[i] = 0, sz[i] = 1;
    }
    
    int root(int x) {
        if (par[x] == x) return x;
        else return par[x] = root(par[x]);
    }
    
    bool issame(int x, int y) {
        return root(x) == root(y);
    }
    
    int size(int x) {
        return sz[root(x)];
    }
    
    bool merge(int x, int y) {
        x = root(x); y = root(y);
        if (x == y) return false;
        if (rank[x] < rank[y]) swap(x, y);
        if (rank[x] == rank[y]) ++rank[x];
        par[y] = x;
        sz[x] += sz[y];
        return true;
    }
};


/* debug */
/*
ostream& operator << (ostream& s, UnionFind uf) {
    map<int, vector<int> > res;
    for (int i = 0; i < uf.par.size(); ++i) {
        int r = uf.root(i);
        res[r].push_back(i);
    }
    for (map<int, vector<int> >::iterator it = res.begin(); it != res.end(); ++it) {
        s << endl;
        for (int j = 0; j < (int)it->second.size(); ++j) {
            s << it->second[j] << ", ";
        }
    }
    return s << endl;
}
 */


int main() {
    int N, Q; cin >> N >> Q;
    UnionFind uf(N);
    for (int q = 0; q < Q; ++q) {
        int p, a, b; cin >> p >> a >> b;
        if (p == 0) uf.merge(a, b);
        else {
            if (uf.issame(a, b)) cout << "Yes" << endl;
            else cout << "No" << endl;
        }
        //cout << uf << endl; // debug
    }
}
