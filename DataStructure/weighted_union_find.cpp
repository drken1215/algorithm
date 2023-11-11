//
// 重み付き Union-Find
//
// reference:
//   drken: 重み付き Union-Find とそれが使える問題のまとめ、および、牛ゲーについて
//     https://qiita.com/drken/items/cce6fc5c579051e64fab
//
// verified:
//   AtCoeer ABC 328 F - Good Set Query
//     https://atcoder.jp/contests/abc328/tasks/abc328_f
//
//   AtCoeer ABC 087 D - People on a Line
//     https://atcoder.jp/contests/abc087/tasks/arc090_b
//


/*
    内部的な変数値 v[0], v[1], ..., v[N-1] を管理する (計算量は基本的に O(α(N)))
    Union-Find 内の各根付き木の根 r に対して v[r] = 0 とする
 
    ・same(x, y): x と y が同じグループにいるかどうか
    ・size(x): x を含むグループの所属メンバー数
    ・groups(): グループ分けの構造を返す, 計算量は O(n)
    ・merge(x, y, w): v[y] - v[x] = w を満たすようにする (すでに同じグループである場合には何もしない)
    ・get_weight(x): v[x] の値を返す
    ・diff(x, y): v[y] - v[x] の値を返す
*/


#include <bits/stdc++.h>
using namespace std;


// Weighted Union-Find (T: the type of v[0], v[1], ..., v[N-1])
template<class T> struct WeightedUnionFind {
    // core member
    vector<int> par;
    vector<T> weight;

    // constructor
    WeightedUnionFind() { }
    WeightedUnionFind(int N, T zero = 0) : par(N, -1), weight(N, zero) {}
    void init(int N, T zero = 0) {
        par.assign(N, -1);
        weight.assign(N, zero);
    }
    
    // core methods
    int root(int x) {
        if (par[x] < 0) return x;
        else {
            int r = root(par[x]);
            weight[x] += weight[par[x]];
            return par[x] = r;
        }
    }
    bool same(int x, int y) {
        return root(x) == root(y);
    }
    int size(int x) {
        return -par[root(x)];
    }
    
    // v[y] - v[x] = w
    bool merge(int x, int y, T w) {
        w += get_weight(x), w -= get_weight(y);
        x = root(x), y = root(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y), w = -w; // merge technique
        par[x] += par[y];
        par[y] = x;
        weight[y] = w;
        return true;
    }
    
    // get v[x]
    T get_weight(int x) {
        root(x);
        return weight[x];
    }
    
    // get v[y] - v[x]
    T get_diff(int x, int y) {
        return get_weight(y) - get_weight(x);
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
    friend ostream& operator << (ostream &s, WeightedUnionFind uf) {
        const vector<vector<int>> &gs = uf.groups();
        for (const vector<int> &g : gs) {
            s << "group: ";
            for (int v : g) s << v << "(" << uf.get_weight(v) << ") ";
            s << endl;
        }
        return s;
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// ABC 328 F
void ABC_328_F() {
    int N, Q;
    cin >> N >> Q;
    WeightedUnionFind<long long> uf(N);
    for (int i = 0; i < Q; ++i) {
        int a, b, d;
        cin >> a >> b >> d;
        --a, --b;
        
        bool good = true;
        if (!uf.same(a, b)) {
            // x[a] - x[b] = d となるように
            uf.merge(b, a, d);
        } else {
            // x[a] - x[b] = d でないとき、ダメ
            if (uf.get_diff(b, a) != d) good = false;
        }
        
        if (good) cout << i+1 << " ";
    }
    cout << endl;
}


// ABC 087 D
void ABC_087_D() {
    int N, M;
    cin >> N >> M;
    WeightedUnionFind<long long> uf(N);
    bool res = true;
    for (int i = 0; i < M; ++i) {
        int l, r, d;
        cin >> l >> r >> d;
        --l, --r;
        if (uf.same(l, r)) {
            if (uf.get_diff(l, r) != d) res = false;
        }
        else uf.merge(l, r, d);
    }
    cout << (res ? "Yes" : "No") << endl;
}


int main() {
    ABC_328_F();
    //ABC_087_D();
}

