//
// 抽象化した全方位木 DP (木 DP パートで辺に関する処理も行う場合)
//
// verified:
//   AtCoder ABC 222 F - Expensive Expense
//     https://atcoder.jp/contests/abc222/tasks/abc222_f
//


#include <bits/stdc++.h>
using namespace std;


// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
template<class T = long long> struct Graph {
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    
    Graph(int n = 0) : list(n), reversed_list(n) { }
    void init(int n = 0) {
        list.assign(n, vector<Edge<T>>());
        reversed_list.assign(n, vector<Edge<T>>());
    }
    Graph &operator = (const Graph &g) {
        list = g.list, reversed_list = g.reversed_list;
        return *this;
    }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    const vector<Edge<T>> &get_rev_edges(int i) const { return reversed_list[i]; }
    const size_t size() const { return list.size(); }
    const void clear() { list.clear(); }
    const void resize(int n) { list.resize(n); }
        
    void add_edge(int from, int to, T val = 1) {
        list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }
    
    void add_bidirected_edge(int from, int to, T val = 1) {
        list[from].push_back(Edge(from, to, val));
        list[to].push_back(Edge(to, from, val));
        reversed_list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }

    friend ostream &operator << (ostream &s, const Graph &G) {
        s << endl;
        for (int i = 0; i < G.size(); ++i) {
            s << i << " -> ";
            for (const auto &e : G[i]) s << e.to << " ";
            s << endl;
        }
        return s;
    }
};  

// 辺に重みがある場合の全方位木 DP
/*
    通常の木 DP において、頂点 v を根とする部分根付き木に関する再帰関数 rec(v) について、
 　　　1. res = IDENTITY
 　　　2. 頂点 v の各子頂点 v2 (その辺を e とする) に対して：res = MERGE(res, ADDEDGE(e, rec(v2)))
 　　　3. return ADDNODE(v, res)
 　　というような更新を行うものとする。
 　　このような木 DP を全方位木 DP へと拡張する。
 */
template<class Monoid, class Weight = long long> struct WeightedReRooting {
    using AddEdgeFunc = function<Monoid(Edge<Weight>, Monoid)>;
    using MergeFunc = function<Monoid(Monoid, Monoid)>;
    using AddNodeFunc = function<Monoid(int, Monoid)>;
    
    // core member
    Graph<Weight> G;  // input graph
    Monoid IDENTITY;
    AddEdgeFunc ADDEDGE;
    MergeFunc MERGE;
    AddNodeFunc ADDNODE;
    
    // inner data
    vector<vector<Monoid>> dp;
    vector<unordered_map<int,int>> ids;
    
    // constructor
    WeightedReRooting() {}
    WeightedReRooting(const Graph<Weight> &g,
                      const AddEdgeFunc &addedge, 
                      const MergeFunc &merge, const AddNodeFunc &addnode, 
                      const Monoid &identity) {
        G = g;
        IDENTITY = identity;
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
            int v2 = G[v][i].to;
            ids[v][v2] = i;
            if (v2 == p) continue;
            dp[v][i] = rec(v2, v);
            res = MERGE(res, ADDEDGE(G[v][i], dp[v][i]));
        }
        return ADDNODE(v, res);
    }
    void rerec(int v, int p, Monoid pval) {
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = G[v][i].to;
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
            int v2 = G[v][i].to, ri = (int)G[v].size() - i - 1;
            if (v2 == p) continue;
            Monoid pval2 = MERGE(left[i], right[ri]);
            rerec(v2, v, ADDNODE(v, pval2));
        }
    }
    void build() {
        dp.assign(G.size(), vector<Monoid>());
        ids.assign(G.size(), unordered_map<int,int>());
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
    Monoid get(int v, int w) {
        return dp[v][ids[v][w]];
    }
    
    // dump
    friend constexpr ostream& operator << (ostream &os, const WeightedReRooting &rr) {
        for (int v = 0; v < rr.G.size(); ++v) {
            for (int i = 0; i < rr.G[v].size(); ++i) {
                os << rr.G[v][i] << ": " << rr.dp[v][i] << endl;
            }
        }
        return os;
    }
};


//------------------------------//
// Examples
//------------------------------//

// AtCoder ABC 222 F - Expensive Expense
void ABC_222_F() {
    long long N, INF = 1LL << 62;
    cin >> N;
    Graph<long long> G(N);
    for (int i = 0; i < N - 1; i++) {
        long long a, b, c; cin >> a >> b >> c, a--, b--;
        G.add_bidirected_edge(a, b, c);
    }
    vector<long long> D(N);
    for (int i = 0; i < N; i++) cin >> D[i];

    auto addedge = [&](Edge<long long> e, long long a) -> long long { 
        return a + e.val;
    };
    auto merge = [&](long long a, long long b) -> long long { 
        return max(a, b);
    };
    auto addnode = [&](int v, long long a) -> long long {
        return max(a, D[v]);
    };
    WeightedReRooting<long long, long long> rr(G, addedge, merge, addnode, 0);

    for (int v = 0; v < N; v++) {
        long long res = 0;
        for (auto e : G[v]) res = max(res, rr.get(v, e.to) + e.val);
        cout << res << endl;
    }
}


int main() {
    ABC_222_F();
}