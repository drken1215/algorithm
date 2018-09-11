//
// Union-Find tree
//
// verified:
//   AOJ Course DSL_1_A Set - Disjoint Set: Union Find Tree
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_1_A&lang=jp
//


/*
    ・併合時の工夫: union by size
    par[x]: x が根のときは x を含むグループのサイズ (の -1 倍)、そうでないときは親ノード

    issame(x, y): x と y が同じグループにいるか, 計算量はならし O(α(n))
    merge(x, y): x と y を同じグループにする, 計算量はならし O(α(n))
    size(x): x を含むグループの所属メンバー数
*/


#include <iostream>
#include <vector>
#include <map>
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
        int com, x, y; cin >> com >> x >> y;
        if (com == 0) uf.merge(x, y);
        else {
            if (uf.issame(x, y)) cout << 1 << endl;
            else cout << 0 << endl;
        }
        //cout << uf << endl;
    }
}
