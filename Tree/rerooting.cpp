//
// 抽象化した全方位木 DP
//
// verified:
//   EDPC V - Subtree
//     https://atcoder.jp/contests/dp/tasks/dp_v
//
//   AtCoder ABC 348 E - Minimize Sum of Distances
//     https://atcoder.jp/contests/abc348/tasks/abc348_e
//

/*
    通常の木 DP において、頂点 v を根とする部分根付き木に関する再帰関数 rec(v) について、
 　　　1. res = IDENTITY
 　　　2. 頂点 v の各子頂点 v2 (その辺を e とする) に対して：res = MERGE(res, rec(v2))
 　　　3. return ADDNODE(v, res)
 　　というような更新を行うものとする。
 　　このような木 DP を全方位木 DP へと拡張する。
 */


#include <bits/stdc++.h>
using namespace std;


// re-rooting
template<class Monoid> struct ReRooting {
    using Graph = vector<vector<int>>;
    using MergeFunc = function<Monoid(Monoid, Monoid)>;
    using AddNodeFunc = function<Monoid(int, Monoid)>;
    
    // core member
    Graph G;
    Monoid IDENTITY;
    MergeFunc MERGE;
    AddNodeFunc ADDNODE;
    
    // inner data
    vector<vector<Monoid>> dp;
    
    // constructor
    ReRooting() {}
    ReRooting(const Graph &g, const Monoid &identity,
              const MergeFunc &merge, const AddNodeFunc &addnode) {
        G = g;
        IDENTITY = identity;
        MERGE = merge;
        ADDNODE = addnode;
        build();
    }
    
    // re-looting dp
    Monoid rec(int v, int p) {
        Monoid res = IDENTITY;
        dp[v].assign(G[v].size(), IDENTITY);
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = G[v][i];
            if (v2 == p) continue;
            dp[v][i] = rec(v2, v);
            res = MERGE(res, dp[v][i]);
        }
        return ADDNODE(v, res);
    }
    void rerec(int v, int p, Monoid pval) {
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = G[v][i];
            if (v2 == p) {
                dp[v][i] = pval;
                continue;
            }
        }
        vector<Monoid> left(G[v].size() + 1, IDENTITY);
        vector<Monoid> right(G[v].size() + 1, IDENTITY);
        for (int i = 0; i < G[v].size(); ++i) {
            left[i + 1] = MERGE(left[i], dp[v][i]);
            right[i + 1] = MERGE(right[i], dp[v][(int)G[v].size() - i - 1]);
        }
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = G[v][i];
            if (v2 == p) continue;
            Monoid pval2 = MERGE(left[i], right[(int)G[v].size() - i - 1]);
            rerec(v2, v, ADDNODE(v, pval2));
        }
    }
    void build() {
        dp.assign(G.size(), vector<Monoid>());
        int root = 0, nullparent = -1;
        rec(root, nullparent);
        rerec(root, nullparent, IDENTITY);
    }
    
    // getter
    Monoid get(int v) {
        Monoid res = IDENTITY;
        for (int i = 0; i < G[v].size(); ++i) {
            res = MERGE(res, dp[v][i]);
        }
        return ADDNODE(v, res);
    }
    
    // dump
    friend constexpr ostream& operator << (ostream &os, const ReRooting<Monoid> &rr) {
        for (int v = 0; v < rr.G.size(); ++v) {
            for (int i = 0; i < rr.G[v].size(); ++i) {
                os << v << " -> " << rr.G[v][i] << ": " << rr.dp[v][i] << endl;
            }
        }
        return os;
    }
};



//------------------------------//
// Examples
//------------------------------//

// TDPC V - Subtree
void TDPC_V() {
    int N, M;
    cin >> N >> M;
    
    using Graph = vector<vector<int>>;
    Graph G(N);
    for (int i = 0; i < N - 1; ++i) {
        int x, y;
        cin >> x >> y;
        --x, --y;
        G[x].push_back(y);
        G[y].push_back(x);
    }
    
    using Monoid = long long;
    Monoid identity = 1;
    auto merge = [&](Monoid a, Monoid b) -> Monoid { return a * b % M; };
    auto addnode = [&](int v, Monoid a) -> Monoid { return (a + 1) % M; };
    ReRooting<Monoid> rr(G, identity, merge, addnode);
    
    //cout << rr << endl;
    
    for (int v = 0; v < N; ++v) {
        cout << (rr.get(v) + M - 1) % M << endl;
    }
}

// ABC 348 E - Minimize Sum of Distances
void ABC_348_E() {
    using Graph = vector<vector<int>>;
    int N;
    cin >> N;
    Graph G(N);
    vector<long long> C(N);
    for (int i = 0; i < N - 1; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        G[a].push_back(b);
        G[b].push_back(a);
    }
    for (int i = 0; i < N; ++i) cin >> C[i];
    
    using Monoid = pair<long long, long long>;  // (siz, sum)
    Monoid IDENTITY = Monoid(-1, -1);
    auto MERGE = [&](Monoid a, Monoid b) {
        if (a.first == -1) return b;
        else if (b.first == -1) return a;
        else return Monoid(a.first + b.first, a.second + b.second);
    };
    auto ADDNODE = [&](int v, Monoid a) {
        Monoid res(C[v], C[v]);
        if (a.first == -1) return res;
        res.first += a.first;
        res.second += a.first + a.second;
        return res;
    };
    
    ReRooting<Monoid> rr(G, IDENTITY, MERGE, ADDNODE);
    long long res = 1LL << 62;
    for (int v = 0; v < N; ++v) {
        auto tmp = rr.get(v);
        res = min(res, tmp.second - tmp.first);
    }
    cout << res << endl;
}


int main() {
    //TDPC_V();
    ABC_348_E();
}

