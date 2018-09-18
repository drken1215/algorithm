//
// maximum stable set O(1.381^{n})
//
// cf.
//   wata: 指数時間アルゴリズム入門
//     https://www.slideshare.net/wata_orz/ss-12131479
//
// verified:
//   CODE THANKS FESTIVAL 2017 G - Mixture Drug
//     https://code-thanks-festival-2017-open.contest.atcoder.jp/tasks/code_thanks_festival_2017_g
//


#include <iostream>
#include <vector>
using namespace std;

// Graph is expressed as adjacent-list
typedef vector<vector<int> > Graph;
int ConnectedCase(Graph &G, vector<int> &can);
int GeneralCase(Graph &G, vector<int> &can);

void dfs(Graph &G, int v, vector<int> &seen, vector<int> &comp, vector<int> &can) {
    seen[v] = true;
    for (auto to : G[v]) {
        if (!seen[to] && can[to]) dfs(G, to, seen, comp, can);
    }
    comp[v] = true;
}

int ConnectedCase(Graph &G, vector<int> &can) {
    int pMax = -1, pMin = -1, Max = -1, Min = (int)G.size(), num = 0;
    for (int i = 0; i < (int)G.size(); ++i) {
        if (!can[i]) continue;
        int tnum = 0; ++num;
        for (int j = 0; j < (int)G[i].size(); ++j) if (can[G[i][j]]) ++tnum;
        if (Max < tnum) Max = tnum, pMax = i;
        if (Min > tnum) Min = tnum, pMin = i;
    }
    if (num == 1) return 1;
    if (Max <= 2) {
        if (Min == 1) return (num+1)/2;
        else return num/2;
    }
    int res = 0;
    if (Min < 2) {
        vector<int> ncan((int)G.size(), 0);
        for (int i = 0; i < (int)G.size(); ++i) ncan[i] = can[i];
        ncan[pMin] = false;
        for (int i = 0; i < (int)G[pMin].size(); ++i) ncan[G[pMin][i]] = false;
        res = max(res, GeneralCase(G, ncan) + 1);
    }
    else {
        vector<int> ncan((int)G.size(), 0);
        for (int i = 0; i < (int)G.size(); ++i) ncan[i] = can[i];
        ncan[pMax] = false;
        for (int i = 0; i < (int)G[pMax].size(); ++i) ncan[G[pMax][i]] = false;
        int temp1 = GeneralCase(G, ncan);
        res = max(res, temp1 + 1);
        
        for (int i = 0; i < (int)G.size(); ++i) ncan[i] = can[i];
        ncan[pMax] = false;
        res = max(res, GeneralCase(G, ncan));
    }
    return res;
}

int GeneralCase(Graph &G, vector<int> &can) {
    if ((int)G.size() == 1) return 1;
    vector<int> seen((int)G.size(), 0);
    int res = 0;
    for (int i = 0; i < (int)G.size(); ++i) {
        if (!seen[i] && can[i]) {
            vector<int> gcan((int)G.size(), false);
            dfs(G, i, seen, gcan, can);
            res += ConnectedCase(G, gcan);
        }
    }
    return res;
}

int StableSet(Graph &G) {
    vector<int> can((int)G.size(), 1);
    return GeneralCase(G, can);
}



int main () {
    int N, M;
    cin >> N >> M;
    Graph G(N);
    for (int i = 0; i < M; ++i) {
        int a, b; cin >> a >> b; --a, --b;
        G[a].push_back(b); G[b].push_back(a);
    }
    cout << StableSet(G) << endl;
}
