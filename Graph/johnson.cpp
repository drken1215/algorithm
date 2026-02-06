//
// Johson のアルゴリズム（全頂点対間最短路を求める, O(EV log V))
//   最初に各頂点のポテンシャル p を求める
//   各辺 e = (u, v) の重みを w(e) + p(u) - p(v) と変更すると、これは非負である
//
// verified
//   AOJ Course GRL_1_C All Pairs Shortest Path
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_1_C&lang=jp
//


#include <bits/stdc++.h>
using namespace std;


template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }

vector<vector<long long>> johnson(const vector<vector<pair<int,long long>>> &G) {
    const long long INF = 1LL << 60;
    int V = G.size();
    vector<vector<long long>> res(V, vector<long long>(V, INF));

    // ポテンシャルを求める
    vector<long long> pot(V, 0);
    for (int iter = 0; iter <= V+1; iter++) {
        bool update = false;
        for (int v = 0; v < V; v++) {
            for (auto [v2, w] : G[v]) {
                if (chmin(pot[v2], pot[v] + w)) update = true;
            }
        }
        if (!update) break;
        if (iter == V && update) return res;  // 負閉路あり
    }

    // 各頂点 s を始点とした Dijkstra 法
    for (int s = 0; s < V; s++) {     
        priority_queue<pair<long long,int>, vector<pair<long long,int>>, greater<pair<long long,int>>> que;
        res[s][s] = 0;
        que.emplace(0, s);
        while (!que.empty()) {
            auto [cur, v] = que.top();
            que.pop();
            if (cur > res[s][v]) continue;
            for (auto [v2, w] : G[v]) {
                long long nw = w + pot[v] - pot[v2];
                assert(nw >= 0);
                if (chmin(res[s][v2], res[s][v] + nw)) que.emplace(res[s][v2], v2);
            }
        }
    }
    
    // ポテンシャル分を補正する
    for (int v = 0; v < V; v++) for (int v2 = 0; v2 < V; v2++) {
        res[v][v2] -= pot[v] - pot[v2];
    }
    return res;
}


//------------------------------//
// Solver
//------------------------------//

// AOJ Course GRL_1_C All Pairs Shortest Path
void AOJ_GRL_1_C() {
    int V, E;
    cin >> V >> E;
    vector<vector<pair<int,long long>>> G(V);
    for (int i = 0; i < E; i++) {
        long long u, v, w;
        cin >> u >> v >> w;
        G[u].emplace_back(v, w);
    }
    const long long INF = 1LL << 60;
    auto res = johnson(G);
    if (res[0][0] != 0) {
        cout << "NEGATIVE CYCLE" << endl;
        return;
    }
    for (int v = 0; v < V; v++) {
        for (int v2 = 0; v2 < V; v2++) {
            if (v2) cout << " ";
            long long dis = res[v][v2];
            if (dis < INF/2) cout << dis;
            else cout << "INF";
        }
        cout << endl;
    }
}


int main() {
    AOJ_GRL_1_C();
}

