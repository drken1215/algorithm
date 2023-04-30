//
// サイクル検出 (by BFS)
//
// cf.
//   BFS (幅優先探索) 超入門！ 〜 キューを鮮やかに使いこなす 〜
//     https://qiita.com/drken/items/996d80bcae64649a6580
//
//
// Verified:
//   AOJ 2891 - な◯りカット
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2891
//


// AOJ 2891 の解答例
#include <iostream>
#include <vector>
#include <queue>
using namespace std;
using Graph = vector<vector<int>>;

int main() {
    // 頂点数 (サイクルを一つ含むグラフなので辺数は N で確定)
    int N;
    cin >> N;

    // グラフ入力受取
    Graph G(N);
    vector<int> deg(N, 0); // 各頂点の出次数
    for (int i = 0; i < N; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b; // 0-indexed にする
        G[a].push_back(b);
        G[b].push_back(a);
        ++deg[a], ++deg[b];
    }

    // 葉たちをキューに挿入する
    // 今回は無向グラフなので、deg[v] == 1 が葉の条件になります
    queue<int> que;
    for (int v = 0; v < N; ++v) if (deg[v] == 1) que.push(v);

    // 探索
    vector<bool> ispushed(N, false);
    while (!que.empty()) {
        int v = que.front(); que.pop();
        ispushed[v] = true;
        for (auto nv : G[v]) {
            --deg[nv];
            if (deg[nv] == 1) que.push(nv);
        }
    }

    // クエリに答える
    int Q; cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        int a, b; cin >> a >> b; --a, --b;
        if (!ispushed[a] && !ispushed[b]) cout << 2 << endl;
        else cout << 1 << endl;
    }
}
