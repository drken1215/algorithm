//
// Weighted Union Find Tree (weight: F2)
//
// cf.
//   https://qiita.com/drken/items/cce6fc5c579051e64fab
//
// verified:
//   Codeforces Round #616 (Div. 1) C. Prefix Enlightenment
//     https://codeforces.com/contest/1290/problem/C
//

/*
  val[v] := 0 or 1 or -1 (-1: non-fixed)
  onum[v] := # of 1 in the connected component
 */


#include <iostream>
#include <vector>
#include <map>
using namespace std;


struct UnionFind {
    vector<int> par;
    vector<int> diff_weight;
    vector<int> val;
    vector<int> onum; // 根が 0 のときの 1 の個数

    UnionFind() { }
    UnionFind(int n) : par(n, -1), diff_weight(n, 0)
                     , val(n, -1), onum(n, 0) {}
    int root(int x) {
        if (par[x] < 0) return x;
        else {
            int r = root(par[x]);
            diff_weight[x] ^= diff_weight[par[x]];
            return par[x] = r;
        }
    }
    
    int calc_weight(int x) {
        int rx = root(x);
        return diff_weight[x];
    }
    
    bool issame(int x, int y) {
        return root(x) == root(y);
    }

    void set(int x, int w) {
        auto rx = root(x);
        auto dw = diff_weight[x];
        val[rx] = w ^ dw;
    }
    
    bool merge(int x, int y, int w = 0) {
        w ^= calc_weight(x); w ^= calc_weight(y);
        x = root(x), y = root(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y); // merge technique
        
        if (w == 0) onum[x] += onum[y];
        else onum[x] += -par[y] - onum[y];

        if (val[y] != -1) val[x] = val[y] ^ w;
        par[x] += par[y];
        par[y] = x;
        diff_weight[y] = w;
        return true;
    }
    
    int diff(int x, int y) {
        return calc_weight(y) ^ calc_weight(x);
    }
    
    int size(int x) {
        return -par[root(x)];
    }

    int get_onum(int x) {
        x = root(x);
        if (val[x] == -1) return min(onum[x], -par[x] - onum[x]);
        else if (val[x] == 0) return onum[x];
        else return -par[x] - onum[x];
    }
};

/* debug */
ostream& operator << (ostream& s, UnionFind uf) {
    map<int, vector<int> > res;
    cout << endl;
    for (int i = 0; i < uf.par.size(); ++i) {
        int r = uf.root(i);
        res[r].push_back(i);
        if (r == i) {
             cout << "root " << r << ": "
                  << uf.val[r] << ", " << uf.onum[r] << endl;
        }
    }
    int iter = 0;
    for (auto it : res) {
        s << (iter++) << "-th group: ";
        for (int j = 0; j < (int)it.second.size(); ++j) {
            s << it.second[j] << "(" << uf.calc_weight(it.second[j]) << "), ";
        }
        s << endl;
    }
    return s << endl;
}


// 入力
int N, K;
string S;
vector<vector<int>> cls;

void solve() {
    long long res = 0;
    UnionFind uf(K);
    
    for (int i = 0; i < N; ++i) {
        int w = 1 - (int)(S[i] - '0');
        int pre_one = -1;
        int cur_one = -1;
        
        if (cls[i].size() == 1) {
            int x = cls[i][0];
            pre_one = uf.get_onum(x);
            uf.set(x, w);
            cur_one = uf.get_onum(x);
            res += (cur_one - pre_one);
        }
        else if (cls[i].size() == 2) {
            int x = cls[i][0], y = cls[i][1];
            if (!uf.issame(x, y)) {
                pre_one = uf.get_onum(x) + uf.get_onum(y);
                uf.merge(x, y, w);
                cur_one = uf.get_onum(x);
                res += (cur_one - pre_one);
            }
        }
        cout << res << endl;
    }
}

int main() {
    cin >> N >> K >> S;
    cls.assign(N, vector<int>());
    for (int k = 0; k < K; ++k) {
        int c; cin >> c;
        for (int i = 0; i < c; ++i) {
            int a; cin >> a; --a;
            cls[a].push_back(k);
        }
    }
    solve();
}
