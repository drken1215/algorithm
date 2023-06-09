//
// Weighted Union Find
//
// cf.
//   https://qiita.com/drken/items/cce6fc5c579051e64fab
//
// verified:
//   ABC 087 D - People on a Line
//     https://beta.atcoder.jp/contests/abc087/tasks/arc090_b
//


#include <bits/stdc++.h>
using namespace std;


// Weighted Union-Find
template<class Abel> struct UnionFind {
    // to be set
    Abel IDENTITY_SUM;
    
    // core member
    vector<int> par;
    vector<Abel> diff_weight;

    // constructor
    UnionFind() { }
    UnionFind(int n, Abel identity = 0) : par(n, -1), diff_weight(n, identity) {}
    void init(int n, Abel identity = 0) {
        par.assign(n, -1);
        diff_weight.assign(n, identity);
    }
    
    // core methods
    int root(int x) {
        if (par[x] < 0) return x;
        else {
            int r = root(par[x]);
            diff_weight[x] += diff_weight[par[x]];
            return par[x] = r;
        }
    }
    bool same(int x, int y) {
        return root(x) == root(y);
    }
    int size(int x) {
        return -par[root(x)];
    }
    
    // calc weight
    Abel calc_weight(int x) {
        root(x);
        return diff_weight[x];
    }
    
    // w[y] - w[x] = w
    bool merge(int x, int y, Abel w) {
        w += calc_weight(x), w -= calc_weight(y);
        x = root(x), y = root(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y), w = -w; // merge technique
        par[x] += par[y];
        par[y] = x;
        diff_weight[y] = w;
        return true;
    }
    
    // calc w[y] - w[x]
    Abel diff(int x, int y) {
        return calc_weight(y) - calc_weight(x);
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
            for (auto v : it.second) {
                s << v << "(" << uf.calc_weight(v) << ") ";
            }
            s << endl;
        }
        return s;
    }
};


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void ABC_087_D() {
    int N, M;
    cin >> N >> M;
    UnionFind<long long> uf(N);
    bool res = true;
    for (int i = 0; i < M; ++i) {
        int l, r, d;
        cin >> l >> r >> d;
        --l, --r;
        if (uf.same(l, r)) {
            long long diff = uf.diff(l, r);
            if (diff != d) res = false;
        }
        else uf.merge(l, r, d);
    }
    cout << (res ? "Yes" : "No") << endl;
}


int main() {
    ABC_087_D();
}






