//
// Partially Persistent Union-Find tree
//
// verified:
//   CODE THANKS FESTIVAL 2017 H - Union Sets
//     https://beta.atcoder.jp/contests/code-thanks-festival-2017-open/tasks/code_thanks_festival_2017_h
//

/*
    par[x]: x が根のときは x を含むグループのサイズ (の -1 倍)、そうでないときは親ノード
    last[x] : x が最後に「根」ではなくなった瞬間の時刻
    history[x] := x の親ノードが変わるたびに (その時刻, 新たな親) を push していく
*/


#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
using namespace std;


using pint = pair<int,int>;
struct PartiallyPersistentUnionFind {
    vector<int> par, last;
    vector<vector<pint> > history;
    
    PartiallyPersistentUnionFind(int n) : par(n, -1), last(n, -1), history(n) {
        for (auto &vec : history) vec.emplace_back(-1, -1);
    }
    void init(int n) {
        par.assign(n, -1); last.assign(n, -1); history.assign(n, vector<pint>());
        for (auto &vec : history) vec.emplace_back(-1, -1);
    }
    
    int root(int t, int x) {
        if (last[x] == -1 || t < last[x]) return x;
        return root(t, par[x]);
    }
    
    bool issame(int t, int x, int y) {
        return root(t, x) == root(t, y);
    }

    bool merge(int t, int x, int y) {
        x = root(t, x); y = root(t, y);
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
};


int main() {
    int N, M, Q; scanf("%d %d", &N, &M);
    PartiallyPersistentUnionFind ppuf(N);
    for (int t = 0; t < M; ++t) {
        int a, b; scanf("%d %d", &a, &b); --a, --b;
        ppuf.merge(t+1, a, b);
    }
    scanf("%d", &Q);
    for (int q = 0; q < Q; ++q) {
        int x, y; scanf("%d %d", &x, &y); --x, --y;
        if (!ppuf.issame(M+10, x, y)) cout << -1 << endl;
        else {
            int low = 0, high = M+10;
            while (high - low > 1) {
                int mid = (low + high) / 2;
                if (ppuf.issame(mid, x, y)) high = mid;
                else low = mid;
            }
            cout << high << endl;
        }
    }
}
