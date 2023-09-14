//
// BFS を用いて、無向グラフの連結成分の個数を求める
//
// cf.
//   BFS (幅優先探索) 超入門！ 〜 キューを鮮やかに使いこなす 〜
//     https://qiita.com/drken/items/996d80bcae64649a6580
//
//

#include <iostream>
#include <vector>
#include <queue>
using namespace std;
using Graph = vector<vector<int>>;

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

    // 頂点 s をスタートとした探索
    vector<int> dist(N, -1);
    queue<int> que;
    int count = 0;
    for (int v = 0; v < N; ++v) {
        if (dist[v] != -1) continue; // v が探索済みならスルー
        dist[v] = 0, que.push(v);
        while (!que.empty()) {
            int v = que.front(); que.pop();
            for (auto nv : G[v]) {
                if (dist[nv] == -1) {
                    dist[nv] = dist[v] + 1;
                    que.push(nv);
                }
            }
        }
        ++count;
    }
    cout << count << endl;
}
