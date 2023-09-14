//
// ツリーの直径を求める (DFS 2 回行う)
//
// verified:
//   AGC 033 C - Removing Coins
//     https://atcoder.jp/contests/agc033/tasks/agc033_c
//


#include <iostream>
#include <vector>
using namespace std;


using Graph = vector<vector<int> >;
struct Diameter {
    vector<int> prev;
    pair<int,int> DiameterDFS(const Graph &G, int v, int p) {
        pair<int,int> res(v, 0);
        for (int i = 0; i < (int)G[v].size(); ++i) {
            if (G[v][i] == p) continue;
            pair<int,int> tmp = DiameterDFS(G, G[v][i], v);
            tmp.second++;
            if (tmp.second > res.second) res = tmp, prev[G[v][i]] = v;
        }
        return res;
    }

    vector<int> solve(const vector<vector<int> > &G) {
        prev.assign((int)G.size(), -1);
        pair<int,int> leaf = DiameterDFS(G, 0, -1);
        prev.assign((int)G.size(), -1);
        pair<int,int> t = DiameterDFS(G, leaf.first, -1);
        vector<int> res;
        int cur = t.first;
        while (cur != -1) res.push_back(cur), cur = prev[cur];
        return res;
    }
};


using Graph = vector<vector<int> >;
int N;
Graph G;

int main() {
    cin >> N;
    G.assign(N, vector<int>());
    for (int i = 0; i < N-1; ++i) {
        int a, b; cin >> a >> b; --a, --b;
        G[a].push_back(b);
        G[b].push_back(a);
    }
    Diameter di;
    auto res = di.solve(G);    
    int V = res.size();
    if (V % 3 == 2) cout << "Second" << endl;
    else cout << "First" << endl;
}
