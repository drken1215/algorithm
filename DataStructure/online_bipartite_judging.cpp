//
// オンライン二部グラフ判定 by 頂点倍加 Union-Find
//
// verified:
//   AtCoder ABC 451 F - Make Bipartite 3
//     https://atcoder.jp/contests/abc451/tasks/abc451_f
//


#include <bits/stdc++.h>
using namespace std;


// オンライン二部グラフ判定 by Union-Find
// left(0 ~ N-1): white, right(N ~ 2n-1): black
struct BipartiteJudge {
    static const int WHITE = 0, BLACK = 1;

    // core member
    int N;
    bool is_bipartite;
    vector<int> par, nex, white, black;

    // constructor
    BipartiteJudge() { }
    BipartiteJudge(int n) : N(n), is_bipartite(true) {
        init(n);
    }
    void init(int n) {
        N = n;
        is_bipartite = true;
        par.assign(n * 2, -1), white.assign(n * 2, 0), black.assign(n * 2, 0);
        nex.resize(n * 2);
        for (int i = 0; i < n * 2; ++i) {
            if (i < n) white[i] = 1;
            else black[i] = 1;
            nex[i] = i;
        }
    }

    // core methods
    bool is_connected(int x, int y) {
        assert(x >= 0 && x < N);
        assert(y >= 0 && y < N);
        return inner_same(x, y) || inner_same(x, y + N);
    }
    bool is_connect_ok(int x, int y) {
        assert(x >= 0 && x < N);
        assert(y >= 0 && y < N);
        return !inner_same(x, y);
    }
    void connect(int x, int y) {
        assert(x >= 0 && x < N);
        assert(y >= 0 && y < N);
        if (is_connected(x, y)) {
            if (!is_connect_ok(x, y)) is_bipartite = false;
            return;
        }
        inner_merge(x, y + N);
        inner_merge(x + N, y);
    }

    // (#white, #black) in case that x is x_color
    pair<int,int> get_nums(int x, bool x_color = WHITE) {
        assert(x >= 0 && x < N);
        if (x_color == BLACK) x += N;
        int r = root(x);
        return make_pair(white[r], black[r]);
    }
    
    // inner methods
    int root(int x) {
        if (par[x] < 0) return x;
        else return par[x] = root(par[x]);
    }
    bool inner_same(int x, int y) {
        return root(x) == root(y);
    }
    bool inner_merge(int x, int y, bool merge_technique = false) {
        x = root(x), y = root(y);
        if (x == y) return false;
        if (merge_technique) if (par[x] > par[y]) swap(x, y); // merge technique
        par[x] += par[y];
        par[y] = x;
        white[x] += white[y], black[x] += black[y];
        swap(nex[x], nex[y]);
        return true;
    }
    
    // get group
    vector<pair<int,int>> group(int x) {
        vector<int> res({x});
        while (nex[res.back()] != x) res.push_back(nex[res.back()]);
        vector<pair<int,int>> ans;
        for (auto &v : res) {
            if (v < N) ans.emplace_back(v, 0);
            else ans.emplace_back(v - N, 1);
        }
        return ans;
    }
    vector<vector<pair<int,int>>> groups() {
        vector<vector<pair<int,int>>> res;
        for (int v = 0; v < N; v++) {
            if (root(v) == v) res.emplace_back(group(v));
        }
        return res;
    }
    
    // debug
    friend ostream& operator << (ostream &s, BipartiteJudge bj) {
        const auto &gs = bj.groups();
        for (const auto &g : gs) {
            s << "group: ";
            for (auto [v, p] : g) s << v << "(" << p << ") ";
            s << endl;
        }
        return s;
    }
};



//------------------------------//
// Examples
//------------------------------//

// AtCoder ABC 451 F - Make Bipartite 3
void ABC_451_F() {
    int N, Q;
    cin >> N >> Q;
    BipartiteJudge bj(N);
    bool ok = true;
    int res = 0;
    for (int q = 0; q < Q; q++) {
        int u, v;
        cin >> u >> v, u--, v--;
        if (bj.is_connected(u, v)) {
            if (!bj.is_connect_ok(u, v)) ok = false;
        } else {
            auto [b1, w1] = bj.get_nums(u);
            auto [b2, w2] = bj.get_nums(v);
            int before = min(b1, w1) + min(b2, w2);
            bj.connect(u, v);
            auto [b, w] = bj.get_nums(u);
            int after = min(b, w);
            res += after - before;
        }
        if (!ok) res = -1;
        cout << res << endl;
    }
}


int main() {
    ABC_451_F();
}