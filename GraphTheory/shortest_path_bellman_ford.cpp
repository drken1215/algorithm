//
// 単一始点最短路 (Bellman-Ford 法, 負辺あり, O(VE))
//
// verified:
//   AOJ Course - Single Source Shortest Path (Negative Edges)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_1_B&lang=ja
//
//   AtCoder ABC 061 D - Score Attack
//     https://atcoder.jp/contests/abc061/tasks/abc061_d
//


#include <bits/stdc++.h>
using namespace std;


// Edge
template<class T> struct Edge {
    int to;
    T weight;
    Edge(int t, T w) : to(t), weight(w) {}
};

// Graph
template<class T> struct Graph : vector<vector<Edge<T>>> {
    Graph() {}
    Graph(int N) : vector<vector<Edge<T>>>(N) {}
};

// Bellman-Ford
template<class T> struct BellmanFord {
    // input
    Graph<T> G;
    
    // results
    // unbounded[v] := whether the length of s-v path is unbounded
    // dp[v] := length of the shortest s-v path
    // prev[v] := previous node of node v on the shortest path
    bool exist_negative = false;
    vector<bool> unbounded;
    vector<T> dp;
    vector<int> prev;
    
    // constructor
    BellmanFord() {}
    BellmanFord(const Graph<T> &graph, int S) { init(graph, S); }
    void init(const Graph<T> &graph, int S) {
        G = graph;
        solve(S);
    }
    
    // find the shortest path from node S
    void solve(int S) {
        // init
        exist_negative = false;
        unbounded.assign(G.size(), false);
        dp.assign(G.size(), numeric_limits<T>::max());
        prev.assign(G.size(), -1);
        
        // bellman-Ford
        dp[S] = 0;
        for (int iter = 0; iter <= G.size()*2; ++iter) {
            for (int v = 0; v < G.size(); ++v) {
                if (dp[v] == numeric_limits<T>::max()) continue;
                for (const auto &e : G[v]) {
                    if (dp[e.to] > dp[v] + e.weight) {
                        dp[e.to] = dp[v] + e.weight;
                        prev[e.to] = v;
                        if (iter == G.size()*2) {
                            exist_negative = true;
                            unbounded[e.to] = true;
                        }
                    }
                }
            }
        }
    }
    
    // reconstruct the shortest path by using prev
    vector<int> reconstruct(int s, int t) {
        vector<int> res;
        int v = t;
        do {
            res.push_back(v);
            v = prev[v];
        } while (v != s && v != -1);
        res.push_back(s);
        reverse(res.begin(), res.end());
        return res;
    }
};


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void AOJ_Single_Source_Shortest_Path_Negative() {
    const long long INF = 1LL<<60;
    int N, M, s;
    cin >> N >> M >> s;
    Graph<long long> G(N);
    for (int i = 0; i < M; ++i) {
        int a, b, c;
        cin >> a >> b >> c;
        G[a].emplace_back(b, c);
    }
    
    BellmanFord<long long> bf(G, s);
    
    if (bf.exist_negative) cout << "NEGATIVE CYCLE" << endl;
    else {
        for (int v = 0; v < N; ++v) {
            if (bf.dp[v] < INF) cout << bf.dp[v] << endl;
            else cout << "INF" << endl;
        }
    }
}

void ABC_061_D() {
    int N, M;
    cin >> N >> M;
    Graph<long long> G(N);
    for (int i = 0; i < M; ++i) {
        int a, b, c;
        cin >> a >> b >> c;
        --a, --b;
        G[a].emplace_back(b, -c);
    }
    
    BellmanFord<long long> bf(G, 0);
    
    // ある負閉路 C が存在して、頂点 0 から C へ到達可能、かつ、C から頂点 N-1 へ到達可能
    if (bf.unbounded[N-1]) cout << "inf" << endl;
    else cout << -bf.dp[N-1] << endl;
}


int main() {
    //AOJ_Single_Source_Shortest_Path_Negative();
    ABC_061_D();
}

