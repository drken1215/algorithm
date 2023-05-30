//
// Union-Find (union by rank)
//
// verified:
//   Yosupo Library Checker - Unionfind
//     https://judge.yosupo.jp/problem/unionfind
//
//   AOJ Course DSL_1_A Set - Disjoint Set: Union Find Tree
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_1_A&lang=jp
//

/*
    併合時の工夫: union by rank
 
    same(x, y): x と y が同じグループにいるか, 計算量はならし O(α(n))
    merge(x, y): x と y を同じグループにする, 計算量はならし O(α(n))
*/


#include <bits/stdc++.h>
using namespace std;


// Union-Find
struct UnionFind {
    // core member
    vector<int> par, rank;

    // constructor
    UnionFind() { }
    UnionFind(int n) : par(n , -1), rank(n, 0) { }
    void init(int n) { par.assign(n, -1), rank.assign(n, 0); }
    
    // core methods
    int root(int x) {
        if (par[x] == -1) return x;
        else return par[x] = root(par[x]);
    }
    
    bool same(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y) {
        x = root(x), y = root(y);
        if (x == y) return false;
        if (rank[x] < rank[y]) swap(x, y);
        if (rank[x] == rank[y]) ++rank[x];
        par[y] = x;
        return true;
    }
    
    // debug
    friend ostream& operator << (ostream &s, UnionFind uf) {
        map<int, vector<int>> groups;
        for (int i = 0; i < uf.par.size(); ++i) {
            int r = uf.root(i);
            groups[r].push_back(i);
        }
        for (const auto &it : groups) {
            s << "group: ";
            for (auto v : it.second) s << v << " ";
            s << endl;
        }
        return s;
    }
};


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void YosupoUnionFind() {
    int N, Q;
    cin >> N >> Q;
    UnionFind uf(N);
    for (int q = 0; q < Q; ++q) {
        int type, x, y;
        cin >> type >> x >> y;
        if (type == 0) uf.merge(x, y);
        else {
            if (uf.same(x, y)) cout << 1 << endl;
            else cout << 0 << endl;
        }
        //cout << uf << endl;
    }
}

void AOJ_DSL_1_A() {
    int N, Q;
    cin >> N >> Q;
    UnionFind uf(N);
    for (int q = 0; q < Q; ++q) {
        int com, x, y;
        cin >> com >> x >> y;
        if (com == 0) uf.merge(x, y);
        else {
            if (uf.same(x, y)) cout << 1 << endl;
            else cout << 0 << endl;
        }
        //cout << uf << endl;
    }
}

    
int main() {
    YosupoUnionFind();
    //AOJ_DSL_1_A();
}


