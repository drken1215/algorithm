//
// 二部グラフ判定 (by DFS)
//
// reference:
//   DFS (深さ優先探索) 超入門！ 〜 グラフ・アルゴリズムの世界への入口 〜【後編】
//     https://qiita.com/drken/items/a803d4fc4a727e02f7ba
//
// verified:
//   AtCoder ARC 327 D - Good Tuple Problem
//     https://atcoder.jp/contests/abc327/tasks/abc327_d
//


#include <bits/stdc++.h>
using namespace std;


// 二部グラフ判定 (color は、1: 黒色, -1: 白色, 0: 未確定)
using Graph = vector<vector<int>>;
bool dfs(const Graph &G, int v, int c, vector<int> &color) {
    color[v] = c;
    for (auto v2 : G[v]) {
        if (color[v2] == 0) {
            if (!dfs(G, v2, -c, color)) return false;
        } else if (color[v2] != -c) return false;
    }
    return true;
}

bool is_bipartite(const Graph &G) {
    int N = (int)G.size();
    vector<int> color(N, 0);
    for (int v = 0; v < N; ++v) {
        if (color[v] != 0) continue;
        if (!dfs(G, v, 1, color)) return false;
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

