//
// Union-Find (union by size)
//
// verified:
//   Yosupo Library Checker - Unionfind
//     https://judge.yosupo.jp/problem/unionfind
//
//   AOJ Course DSL_1_A Set - Disjoint Set: Union Find Tree
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_1_A&lang=jp
//


/*
    same(x, y): x と y が同じグループにいるか, 計算量はならし O(α(n))
    merge(x, y): x と y を同じグループにする, 計算量はならし O(α(n))
    size(x): x を含むグループの所属メンバー数, 計算量はならし O(α(n))
    groups(): グループ分けの構造を返す, 計算量は O(n)
*/


#include <bits/stdc++.h>
using namespace std;


// Union-Find
struct UnionFind {
    // core member
    vector<int> par;

    // constructor
    UnionFind() { }
    UnionFind(int n) : par(n, -1) { }
    void init(int n) { par.assign(n, -1); }
    
    // core methods
    int root(int x) {
        if (par[x] < 0) return x;
        else return par[x] = root(par[x]);
    }
    
    bool same(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y) {
        x = root(x), y = root(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y); // merge technique
        par[x] += par[y];
        par[y] = x;
        return true;
    }
    
    int size(int x) {
        return -par[root(x)];
    }
    
    // get groups
    vector<vector<int>> groups() {
        vector<vector<int>> member(par.size());
        for (int v = 0; v < (int)par.size(); ++v) {
            member[root(v)].push_back(v);
        }
        vector<vector<int>> res;
        for (int v = 0; v < (int)par.size(); ++v) {
            if (!member[v].empty()) res.push_back(member[v]);
        }
        return res;
    }
    
    // debug
    friend ostream& operator << (ostream &s, UnionFind uf) {
        const vector<vector<int>> &gs = uf.groups();
        for (const vector<int> &g : gs) {
            s << "group: ";
            for (int v : g) s << v << " ";
            s << endl;
        }
        return s;
    }
};



//------------------------------//
// Examples
//------------------------------//

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
        cout << uf << endl;
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

