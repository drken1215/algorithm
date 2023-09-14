//
// 二部グラフ判定 (by DFS)
//
// cf.
//   DFS (深さ優先探索) 超入門！ 〜 グラフ・アルゴリズムの世界への入口 〜【後編】
//     https://qiita.com/drken/items/a803d4fc4a727e02f7ba
//
//
// verified:
//   ARC 099 E - Independence
//     https://beta.atcoder.jp/contests/arc099/tasks/arc099_c
//

/*
    二部グラフかどうかを判定する
    また、各連結成分ごとに、(左ノード数, 右ノード数) を求める
 */


#include <iostream>
#include <vector>
using namespace std;


// nums[i] := i 番目の連結成分の {左ノード数, 右ノード数}, 「dir = 1: 左、-1: 右」
using pint = pair<int,int>;
using Graph = vector<vector<int> >;

bool dfs(const Graph &G, int v, int vdir, pint &num, vector<int> &dir) {
    bool res = true;
    dir[v] = vdir;
    if (vdir == 1) ++num.first;
    else if (vdir == -1) ++num.second;
    for (auto nv : G[v]) {
        if (dir[nv] == 0) {
            if (!dfs(G, nv, -vdir, num, dir)) res = false;
        }
        else if (dir[nv] != -vdir) res = false;
    }
    return res;
}

bool isbipartite(const Graph &G, vector<pint> &nums) {
    bool res = true;
    int N = (int)G.size();
    vector<int> dir(N, 0);
    for (int v = 0; v < N; ++v) {
        if (dir[v] != 0) continue;
        pint num = {0, 0};
        if (!dfs(G, v, 1, num, dir)) res = false;
        nums.push_back(num);
    }
    return res;
}


int main() {
    int N, M;
    cin >> N >> M;
    vector<vector<int> > isexistedge(N, vector<int>(N, 1));
    for (int i = 0; i < M; ++i) {
        int a, b;
        cin >> a >> b; --a, --b;
        isexistedge[a][b] = isexistedge[b][a] = 0;
    }
    Graph G(N);
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j)
        if (isexistedge[i][j] && i != j)
            G[i].push_back(j);
    
    // 二部グラフの構成
    vector<pint> nums;
    if (!isbipartite(G, nums)) {
        cout << -1 << endl;
        return 0;
    }
    
    // ナップサック DP
    vector<int> dp(N+10, 0);
    dp[0] = 1;
    for (int i = 0; i < (int)nums.size(); ++i) {
        for (int j = N; j >= 0; --j) {
            if (!dp[j]) continue;
            dp[j] = 0;
            dp[j + nums[i].first] = 1;
            dp[j + nums[i].second] = 1;
        }
    }
    
    long long a = -1;
    long long mindif = N;
    for (int i = 0; i <= N; ++i) {
        if (!dp[i]) continue;
        if (mindif > abs(i - N/2)) mindif = abs(i - N/2), a = i;
    }
    long long b = N-a;
    cout << a*(a-1)/2 + b*(b-1)/2 << endl;
}
