//
// 重みが 0 と 1 のみであるグラフの最短路問題 (by 0-1 BFS)
//
// References:
//   AtCoder ARC 005 C - 器物損壊！高橋君 (試験管水色)
//     https://drken1215.hatenablog.com/entry/2021/07/30/024800
//
// verified:
//   AtCoder ARC 005 C - 器物損壊！高橋君
//     https://atcoder.jp/contests/arc005/tasks/arc005_3
//


#include <iostream>
#include <vector>
#include <queue>
#include <string>
using namespace std;

// 無限大を表す値
const int INF = 1<<29;

// 上下左右への動きを表すベクトル
const vector<int> dx = {1, 0, -1, 0};
const vector<int> dy = {0, 1, 0, -1};

int main() {
    // 入力
    int H, W;  // 縦の長さ, 横の長さ
    cin >> H >> W;
    vector<string> field(H);
    for (int i = 0; i < H; ++i) cin >> field[i];

    // スタートとゴールのマスを割り出す
    int sx = -1, sy = -1, gx = -1, gy = -1;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            if (field[i][j] == 's') sx = i, sy = j;
            if (field[i][j] == 'g') gx = i, gy = j;
        }
    }

    // 各頂点は pair<int,int> 型で表すことにする
    using Node = pair<int,int>;
    deque<Node> que;  // deque

    // 初期条件
    // dist[i][j] はマス (i, j) への最短路長を表す
    que.push_front(Node(sx, sy));
    vector<vector<int>> dist(H, vector<int>(W, INF));
    dist[sx][sy] = 0;

    // 0-1 BFS 開始
    while (!que.empty()) {
        // deque の先頭の要素を取り出す
        auto [x, y] = que.front();
        que.pop_front();

        // 隣接頂点を順にみていく
        for (int dir = 0; dir < 4; ++dir) {
            int nx = x + dx[dir];
            int ny = y + dy[dir];

            // 盤面外はスキップ
            if (nx < 0 || nx >= H || ny < 0 || ny >= W) continue;

            // コスト 0 の場合は、deque の先頭に push
            if (field[nx][ny] != '#') {
                if (dist[nx][ny] > dist[x][y]) {
                    dist[nx][ny] = dist[x][y];
                    que.push_front(Node(nx, ny));
                }
            }
            // コスト 1 の場合は、deque の末尾に push
            else {
                if (dist[nx][ny] > dist[x][y] + 1) {
                    dist[nx][ny] = dist[x][y] + 1;
                    que.push_back(Node(nx, ny));
                }
            }
        }
    }

    // 最短路長
    if (dist[gx][gy] <= 2) cout << "YES" << endl;
    else cout << "NO" << endl;
}
