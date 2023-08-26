//
// 単一始点最短路 (Dijkstra 法, 正辺のみ, O(V + E log V))
//
// verified:
//   Yosupo Judge - Shortest Path
//     https://judge.yosupo.jp/problem/shortest_path
//
//   AOJ Course - Single Source Shortest Path
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_1_A&lang=ja
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

// Dijkstra, find {dp, prev}
// dp[v] := length of the shortest s-v path
// prev[v] := previous node of node v on the shortest path
template<class T> pair<vector<T>, vector<int>> Dijkstra(const Graph<T> &G, int s) {
    // Dijkstra 法実行のための変数
    vector<T> dp(G.size(), numeric_limits<T>::max());
    vector<int> prev(G.size(), -1);
    using Node = pair<T, int>;
    priority_queue<Node, vector<Node>, greater<Node>> que;
    
    // 初期条件
    dp[s] = 0;
    que.push(Node(0, s));
    
    // Dijkstra 法
    while (!que.empty()) {
        const auto [cur, v] = que.top();
        que.pop();
        
        // 枝刈り
        if (cur > dp[v]) continue;
        
        // 各辺への遷移を試す
        for (const auto &e : G[v]) {
            if (dp[e.to] > dp[v] + e.weight) {
                dp[e.to] = dp[v] + e.weight;
                prev[e.to] = v;
                que.push(Node(dp[e.to], e.to));
            }
        }
    }
    return make_pair(dp, prev);
}

// reconstruct the shortest path by using prev
vector<int> reconstruct(const vector<int> &prev, int s, int t) {
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


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/
    
void Yosupo_Shortest_Path() {
    int N, M, s, t;
    cin >> N >> M >> s >> t;
    Graph<long long> G(N);
    for (int i = 0; i < M; ++i) {
        int a, b, c;
        cin >> a >> b >> c;
        G[a].emplace_back(b, c);
    }
    
    const auto [dp, prev] = Dijkstra(G, s);
    const auto &path = reconstruct(prev, s, t);
    
    if (prev[t] == -1) cout << -1 << endl;
    else {
        cout << dp[t] << " " << path.size()-1 << endl;
        for (int i = 0; i+1 < path.size(); ++i) {
            cout << path[i] << " " << path[i+1] << endl;
        }
    }
}

void AOJ_Single_Source_Shortest_Path() {
    const long long INF = 1LL<<60;
    int N, M, s;
    cin >> N >> M >> s;
    Graph<long long> G(N);
    for (int i = 0; i < M; ++i) {
        int a, b, c;
        cin >> a >> b >> c;
        G[a].emplace_back(b, c);
    }
    const auto [dp, prev] = Dijkstra(G, s);
    for (int v = 0; v < N; ++v) {
        if (dp[v] < INF) cout << dp[v] << endl;
        else cout << "INF" << endl;
    }
}


int main() {
    Yosupo_Shortest_Path();
    //AOJ_Single_Source_Shortest_Path();
}


