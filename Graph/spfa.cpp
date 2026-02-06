//
// SPFA（shortest path faster algorithm）
//   reference: https://hogloid.hatenablog.com/entry/20120409/1333973448
//
// verified
//   AOJ Course GRL_1_B - Single Source Shortest Path (Negative Edges)
//     https://onlinejudge.u-aizu.ac.jp/courses/library/5/GRL/1/GRL_1_B
//


#include <bits/stdc++.h>
using namespace std;


// 負閉路があるかどうかも判定する
pair<bool, vector<long long>> spfa(const vector<vector<pair<int,long long>>> &G, int s) {
    const long long INF = 1LL << 60;
    int V = (int)G.size();
    vector<long long> res(V, INF);
    queue<int> que;
    vector<bool> inque(V, false);
    vector<int> cnt(V, 0);
    que.push(s), res[s] = 0, inque[s] = true;
    while (!que.empty()) {
        int cur = que.front();
        que.pop();
        inque[cur] = false;
        if (cnt[cur] > V) return {false, res};  // include negative-cycle
        cnt[cur]++;
        for (auto [nex, w]: G[cur]) {
            if (res[nex] > res[cur] + w) {
                res[nex] = res[cur] + w;
                if (!inque[nex]) inque[nex] = true, que.push(nex);
            }
        }
    }
    return {true, res};
}


//------------------------------//
// Solver
//------------------------------//

// AOJ Course GRL_1_B - Single Source Shortest Path (Negative Edges)
void AOJ_GRL_1_B() {
    long long INF = 1LL << 60;
    int V, E, r;
    cin >> V >> E >> r;
    vector<vector<pair<int,long long>>> G(V);
    for (int i = 0; i < E; i++) {
        long long s, t, d;
        cin >> s >> t >> d;
        G[s].emplace_back(t, d);
    }
    auto [flag, res] = spfa(G, r);
    if (flag) {
        for (int v = 0; v < V; v++) {
            if (res[v] < INF/2) cout << res[v] << endl;
            else cout << "INF" << endl;
        }
    } else {
        cout << "NEGATIVE CYCLE" << endl;
    }
}


int main() {
    AOJ_GRL_1_B();
}

