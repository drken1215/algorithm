//
// 抽象化した全方位木 DP
//
// verified:
//   AtCoder ABC 348 E - Minimize Sum of Distances
//     https://atcoder.jp/contests/abc348/tasks/abc348_e
//

/*
    通常の木 DP において、頂点 v を根とする部分根付き木に関する再帰関数 rec(v) について、
 　　　1. res = IDENTITY
 　　　2. 頂点 v の各子頂点 v2 (その辺を e とする) に対して：res = MERGE(res, ADDEDGE(e, rec(v2)))
 　　　3. return ADDNODE(v, res)
 　　というような更新を行うものとする。
 　　このような木 DP を全方位木 DP へと拡張する。
 */


#include <bits/stdc++.h>
using namespace std;


// re-rooting
template<class Monoid, class Edge> struct ReRooting {
    using Graph = vector<vector<Edge>>;
    using GetIdFunc = function<int(Edge)>;
    using AddEdgeFunc = function<Monoid(Edge, Monoid)>;
    using MergeFunc = function<Monoid(Monoid, Monoid)>;
    using AddNodeFunc = function<Monoid(int, Monoid)>;
    
    // core member
    Graph G;
    Monoid IDENTITY;
    GetIdFunc GETID;
    AddEdgeFunc ADDEDGE;
    MergeFunc MERGE;
    AddNodeFunc ADDNODE;
    
    // inner data
    vector<vector<Monoid>> dp;
    
    // constructor
    ReRooting() {}
    ReRooting(const Graph &g, const Monoid &identity, const GetIdFunc &getid,
              const AddEdgeFunc &addedge, const MergeFunc &merge, const AddNodeFunc &addnode) {
        G = g;
        IDENTITY = identity;
        GETID = getid;
        ADDEDGE = addedge;
        MERGE = merge;
        ADDNODE = addnode;
        build();
    }
    
    // re-looting dp
    Monoid rec(int v, int p) {
        Monoid res = IDENTITY;
        dp[v].assign(G[v].size(), IDENTITY);
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = GETID(G[v][i]);
            if (v2 == p) continue;
            dp[v][i] = rec(v2, v);
            res = MERGE(res, ADDEDGE(G[v][i], dp[v][i]));
        }
        return ADDNODE(v, res);
    }
    void rerec(int v, int p, Monoid pval) {
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = GETID(G[v][i]);
            if (v2 == p) {
                dp[v][i] = pval;
                continue;
            }
        }
        vector<Monoid> left(G[v].size() + 1, IDENTITY);
        vector<Monoid> right(G[v].size() + 1, IDENTITY);
        for (int i = 0; i < G[v].size(); ++i) {
            int ri = (int)G[v].size() - i - 1;
            left[i + 1] = MERGE(left[i], ADDEDGE(G[v][i], dp[v][i]));
            right[i + 1] = MERGE(right[i], ADDEDGE(G[v][ri], dp[v][ri]));
        }
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = GETID(G[v][i]), ri = (int)G[v].size() - i - 1;
            if (v2 == p) continue;
            Monoid pval2 = MERGE(left[i], right[ri]);
            rerec(v2, v, ADDNODE(v, pval2));
        }
    }
    void build() {
        dp.assign(G.size(), vector<Monoid>());
        int root = 0;
        rec(root, -1);
        rerec(root, -1, IDENTITY);
    }
    
    // getter
    Monoid get(int v) {
        Monoid res = IDENTITY;
        for (int i = 0; i < G[v].size(); ++i) {
            res = MERGE(res, ADDEDGE(G[v][i], dp[v][i]));
        }
        return ADDNODE(v, res);
    }
};



//------------------------------//
// Examples
//------------------------------//

// ABC 348 E
void ABC_348_E() {
    using Edge = int;
    using Graph = vector<vector<Edge>>;
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
    auto GETID = [&](Edge e) { return e; };
    auto ADDEDGE = [&](Edge e, Monoid a) {
        return a;
    };
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
    
    ReRooting<Monoid, Edge> rr(G, IDENTITY, GETID, ADDEDGE, MERGE, ADDNODE);
    long long res = 1LL << 62;
    for (int v = 0; v < N; ++v) {
        auto tmp = rr.get(v);
        res = min(res, tmp.second - tmp.first);
    }
    cout << res << endl;
}



int main() {
    ABC_348_E();
}

