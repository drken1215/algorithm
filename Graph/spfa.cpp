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
pair<bool, vector<long long> spfa(const vector<vector<pair<int,long long>>> &G) {
    queue<int> que;
    vector<bool> inque(size(), true);
    vector<int> cnt(size(), 0);
    for (int v = 0; v < size(); v++) que.push(v);
    while (!que.empty()) {
        int cur = que.front();
        que.pop();
        inque[cur] = false;
        if (cnt[cur] > size()) return false;  // include negative-cycle
        cnt[cur]++;
        for (const auto &e : G[cur]) {
            if (!e.cap) continue;
            if (pot[e.to] > pot[cur] + e.cost) {
                pot[e.to] = pot[cur] + e.cost;
                if (!inque[e.to]) inque[e.to] = true, que.push(e.to);
            }
        }
    }
    return true;
}


//------------------------------//
// Solver
//------------------------------//

// AOJ Course GRL_1_B - Single Source Shortest Path (Negative Edges)
void AOJ_GRL_1_B() {
    
}


int main() {
    AOJ_GRL_1_B();
}

