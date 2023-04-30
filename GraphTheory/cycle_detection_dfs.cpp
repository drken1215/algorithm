//
// トポロジカルソート (by BFS)
//
// cf.
//   DFS (深さ優先探索) 超入門！ 〜 グラフ・アルゴリズムの世界への入口 〜【後編】
//     https://qiita.com/drken/items/a803d4fc4a727e02f7ba
//
//
// Verified:
//   AOJ 2891 - な◯りカット
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2891
//


// AOJ 2891 の解答例
#include <iostream>
#include <vector>
#include <set>
#include <stack>
using namespace std;
using Graph = vector<vector<int>>;

// 探索
vector<bool> seen, finished;

// サイクル復元のための情報
int pos = -1; // サイクル中に含まれる頂点 pos
stack<int> hist; // 訪問履歴

void dfs(const Graph &G, int v, int p) {
    seen[v] = true;
    hist.push(v);
    for (auto nv : G[v]) {
        if (nv == p) continue; // 逆流を禁止する

        // 完全終了した頂点はスルー
        if (finished[nv]) continue;

        // サイクルを検出
        if (seen[nv] && !finished[nv]) {
            pos = nv;
            return;
        }

        // 再帰的に探索
        dfs(G, nv, v);

        // サイクル検出したならば真っ直ぐに抜けていく
        if (pos != -1) return;
    }
    hist.pop();
    finished[v] = true;
}

int main() {
    // 頂点数 (サイクルを一つ含むグラフなので辺数は N で確定)
    int N; cin >> N;

    // グラフ入力受取
    Graph G(N);
    for (int i = 0; i < N; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b; // 頂点番号が 1-indexed で与えられるので 0-indexed にする
        G[a].push_back(b);
        G[b].push_back(a);
    }

    // 探索
    seen.assign(N, false), finished.assign(N, false);
    pos = -1;
    dfs(G, 0, -1);

    // サイクルを復元
    set<int> cycle;
    while (!hist.empty()) {
        int t = hist.top();
        cycle.insert(t);
        hist.pop();
        if (t == pos) break;
    }

    // クエリに答える
    int Q; cin >> Q;
    for (int _ = 0; _ < Q; ++_) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        if (cycle.count(a) && cycle.count(b)) cout << 2 << endl;
        else cout << 1 << endl;
    }
}
