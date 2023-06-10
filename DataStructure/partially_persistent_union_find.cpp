//
// Partially Persistent Union-Find
//
// verified:
//   AGC 002 D - Stamp Rally
//     https://beta.atcoder.jp/contests/agc002/tasks/agc002_d
//

/*
    t 秒後の状態を参照できる Union-Find 木です。
    t 秒後の状態からの再構築はできなくて、それができるものは「完全永続 Union-Find 木」と呼ばれます。
 
    par[x]: x が根のときは x を含むグループのサイズ (の -1 倍)、そうでないときは親ノード
    last[x] : x が最後に「根」ではなくなった瞬間の時刻
    history[x] := x の親ノードが変わるたびに (その時刻, 新たな親) を push していく
*/


#include <bits/stdc++.h>
using namespace std;


// partially persistent union-find
struct PartiallyPersistentUnionFind {
    // core member
    vector<int> par, last;
    vector<vector<pair<int,int>>> history;

    // constructor
    PartiallyPersistentUnionFind() { }
    PartiallyPersistentUnionFind(int n) : par(n, -1), last(n, -1), history(n) {
        for (auto &vec : history) vec.emplace_back(-1, -1);
    }
    
    // core methods
    void init(int n) {
        par.assign(n, -1);
        last.assign(n, -1);
        history.assign(n, vector<pair<int,int>>());
        for (auto &vec : history) vec.emplace_back(-1, -1);
    }
    
    int root(int t, int x) {
        if (last[x] == -1 || t < last[x]) return x;
        return root(t, par[x]);
    }
    
    bool same(int t, int x, int y) {
        return root(t, x) == root(t, y);
    }
    
    bool merge(int t, int x, int y) {
        x = root(t, x), y = root(t, y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y); // merge technique
        par[x] += par[y];
        par[y] = x;
        last[y] = t;
        history[x].emplace_back(t, par[x]);
        return true;
    }
    
    int size(int t, int x) {
        x = root(t, x);
        return -prev(lower_bound(history[x].begin(), history[x].end(), make_pair(t, 0)))->second;
    }
    
    // debug
    friend ostream& operator << (ostream &s, PartiallyPersistentUnionFind uf) {
        map<int, vector<int>> groups;
        int lasttime = 0;
        for (int i = 0; i < uf.par.size(); ++i) {
            lasttime = max(lasttime, uf.last[i]);
        }
        for (int i = 0; i < uf.par.size(); ++i) {
            int r = uf.root(lasttime+1, i);
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

void AGC_002_D() {
    int N, M, Q;
    cin >> N >> M;
    PartiallyPersistentUnionFind uf(N);
    for (int t = 0; t < M; ++t) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        uf.merge(t, a, b);
    }
    cin >> Q;
    for (int q = 0; q < Q; ++q) {
        int x, y, z;
        cin >> x >> y >> z;
        --x, --y;
        int low = -1, high = M;
        while (high - low > 1) {
            int mid = (low + high) / 2;
            int num = (uf.same(mid, x, y) ? uf.size(mid, x) : uf.size(mid, x) + uf.size(mid, y));
            if (num >= z) high = mid;
            else low = mid;
        }
        cout << high+1 << endl;
    }
}


int main() {
    AGC_002_D();
}


