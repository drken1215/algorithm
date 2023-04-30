//
// DFS を用いて、無向グラフの連結成分の個数を求める
//
// cf.
//   DFS (深さ優先探索) 超入門！ 〜 グラフ・アルゴリズムの世界への入口 〜【後編】
//     https://qiita.com/drken/items/a803d4fc4a727e02f7ba
//

#include <iostream>
#include <vector>
using namespace std;
using Graph = vector<vector<int>>;

// 深さ優先探索
vector<bool> seen;
void dfs(const Graph &G, int v) {
    seen[v] = true;
    for (auto next_v : G[v]) {
        if (seen[next_v]) continue;
        dfs(G, next_v); // 再帰的に探索
    }
}

int main() {
    // 頂点数と辺数
    int N, M; cin >> N >> M;

    // グラフ入力受取
    Graph G(N);
    for (int i = 0; i < M; ++i) {
        int a, b;
        cin >> a >> b;
        G[a].push_back(b);
        G[b].push_back(a);
    }

    // 全頂点が訪問済みになるまで探索
    int count = 0;
    seen.assign(N, false);
    for (int v = 0; v < N; ++v) {
        if (seen[v]) continue; // v が探索済みだったらスルー
        dfs(G, v); // v が未探索なら v を始点とした DFS を行う
        ++count;
    }
    cout << count << endl;
}
