//
// Weighted Union Find Tree
//
// c.f.
//   https://qiita.com/drken/items/cce6fc5c579051e64fab
//
// verified:
//   ABC 087 D - People on a Line
//     https://beta.atcoder.jp/contests/abc087/tasks/arc090_b
//


/*
 ・併合時の工夫: union by size
    par[x]: x が根のときは x を含むグループのサイズ (の -1 倍)、そうでないときは親ノード
 
*/


#include <iostream>
#include <vector>
#include <map>
using namespace std;


template<class Abel> struct UnionFind {
    const Abel UNITY_SUM = 0;      // to be set
    vector<int> par;
    vector<Abel> diff_weight;
    
    UnionFind(int n) : par(n, -1), diff_weight(n, UNITY_SUM) {}
    void init(int n) { par.assign(n, -1), diff_weight.assign(n, UNITY_SUM); }
    
    int root(int x) {
        if (par[x] < 0) return x;
        else {
            int r = root(par[x]);
            diff_weight[x] += diff_weight[par[x]];
            return par[x] = r;
        }
    }
    
    Abel calc_weight(int x) {
        root(x);
        return diff_weight[x];
    }
    
    bool issame(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y, Abel w = 0) {
        w += calc_weight(x); w -= calc_weight(y);
        x = root(x); y = root(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y), w = -w; // merge technique
        par[x] += par[y];
        par[y] = x;
        diff_weight[y] = w;
        return true;
    }
    
    Abel diff(int x, int y) {
        return calc_weight(y) - calc_weight(x);
    }
    
    int size(int x) {
        return -par[root(x)];
    }
};


/* debug */
/*
template<class Abel> ostream& operator << (ostream& s, UnionFind<Abel> uf) {
    map<int, vector<int> > res;
    for (int i = 0; i < uf.par.size(); ++i) {
        int r = uf.root(i);
        res[r].push_back(i);
    }
    for (map<int, vector<int> >::iterator it = res.begin(); it != res.end(); ++it) {
        s << endl;
        for (int j = 0; j < (int)it->second.size(); ++j) {
            s << it->second[j] << "(" << uf.calc_weight(it->second[j]) << "), ";
        }
    }
    return s << endl;
}
*/


int main() {
    int N, M;
    cin >> N >> M;
    UnionFind<long long> uf(N);
    for (int i = 0; i < M; ++i) {
        int l, r, d; cin >> l >> r >> d; --l, --r;
        if (uf.issame(l, r)) {
            long long diff = uf.diff(l, r);
            if (diff != d) {
                puts("No");
                return 0;
            }
        }
        else {
            uf.merge(l, r, d);
        }
    }
    puts("Yes");
}
