//
// 二部グラフ判定 (by BFS)
//
// reference:
//   BFS (幅優先探索) 超入門！ 〜 キューを鮮やかに使いこなす 〜
//     https://qiita.com/drken/items/996d80bcae64649a6580
//
// verified:
//   AtCoder ARC 327 D - Good Tuple Problem
//     https://atcoder.jp/contests/abc327/tasks/abc327_d
//


#include <bits/stdc++.h>
using namespace std;


// 二部グラフ判定 (color は、1: 黒色, -1: 白色, 0: 未確定)
using Graph = vector<vector<int>>;
bool is_bipartite(const Graph &G) {
    int N = (int)G.size();
    vector<int> color(N, 0);
    for (int v = 0; v < N; ++v) {
        if (color[v] != 0) continue;
        color[v] = 1;
        queue<int> que;
        que.push(v);
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            for (auto v2 : G[v]) {
                if (color[v2] != 0) {
                    if (color[v2] == color[v]) return false;
                } else {
                    color[v2] = -color[v];
                    que.push(v2);
                }
            }
        }
    }
    return true;
}



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void ABC_327_D() {
    int N, M;
    cin >> N >> M;
    vector<int> A(M), B(M);
    for (int i = 0; i < M; ++i) cin >> A[i], --A[i];
    for (int i = 0; i < M; ++i) cin >> B[i], --B[i];
    
    vector<vector<int>> G(N);
    for (int i = 0; i < M; ++i) {
        G[A[i]].push_back(B[i]);
        G[B[i]].push_back(A[i]);
    }
    cout << (is_bipartite(G) ? "Yes" : "No") << endl;
}


int main() {
    ABC_327_D();
}

