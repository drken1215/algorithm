//
// Union-Find において、要素 x を含むグループを取得できるようにしたもの
//
// verified:
//   AtCoder ABC 328 F - Good Set Query
//     https://atcoder.jp/contests/abc328/tasks/abc328_f
//


/*
    same(x, y): x と y が同じグループにいるか, 計算量はならし O(α(N))
    merge(x, y): x と y を同じグループにする, 計算量はならし O(α(N))
    size(x): x を含むグループの所属メンバー数, 計算量はならし O(α(N))
    group(x): x を含むグループのメンバーを返す, 計算量はグループのサイズを s として O(s)
    groups(): グループ分けの構造を返す, 計算量は O(N)
*/


#include <bits/stdc++.h>
using namespace std;


// Union-Find
struct UnionFind {
    // core member
    vector<int> par, nex;

    // constructor
    UnionFind() { }
    UnionFind(int N) : par(N, -1), nex(N) {
        init(N);
    }
    void init(int N) {
        par.assign(N, -1);
        nex.resize(N);
        for (int i = 0; i < N; ++i) nex[i] = i;
    }
    
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
        swap(nex[x], nex[y]);
        return true;
    }
    
    int size(int x) {
        return -par[root(x)];
    }
    
    // get group
    vector<int> group(int x) {
        vector<int> res({x});
        while (nex[res.back()] != x) res.push_back(nex[res.back()]);
        return res;
    }
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



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// ABC 328 F
void ABC_328_F() {
    int N, Q;
    cin >> N >> Q;
    UnionFind uf(N);
    vector<long long> x(N, 0);  // x value
    for (int i = 0; i < Q; ++i) {
        int a, b, d;
        cin >> a >> b >> d;
        --a, --b;
        
        bool good = true;
        if (uf.same(a, b)) {
            // x[a] - x[b] = d でないとき、ダメ
            if (x[a] - x[b] != d) good = false;
        } else {
            // マージテクにより、size(a) < size(b) となるようにする
            if (uf.size(a) > uf.size(b)) {
                swap(a, b);
                d = -d;
            }
            
            // a を含むグループのメンバーたちに足す値を求める
            long long add = (x[b] + d) - x[a];
            
            // x の値を調整して、マージする
            vector<int> ids = uf.group(a);
            for (int id : ids) x[id] += add;
            uf.merge(a, b);
        }

        // good ならば出力する
        if (good) cout << i+1 << " ";
    }
    cout << endl;
}


int main() {
    ABC_328_F();
}
