//
// 木の直径を求める
//
// verified:
//   Yosupo Library Checker - Tree Diameter
//     https://judge.yosupo.jp/problem/tree_diameter
//
//   AGC 033 C - Removing Coins
//     https://atcoder.jp/contests/agc033/tasks/agc033_c
//


#include <bits/stdc++.h>
using namespace std;


// find diameter of tree
template<class Graph = vector<vector<int>>> struct Diameter {
    vector<int> path, prev;

    Diameter() {}
    Diameter(const Graph &G) {
        solve(G);
    }
    pair<int, int> DiameterDFS(const Graph &G, int v, int p) {
        pair<int, int> res(v, 0);
        for (auto to : G[v]) {
            if (to == p) continue;
            pair<int, int> tmp = DiameterDFS(G, to, v);
            tmp.second++;
            if (tmp.second > res.second) res = tmp, prev[to] = v;
        }
        return res;
    }
    vector<int> solve(const Graph &G) {
        prev.assign((int)G.size(), -1);
        auto [leaf, distance] = DiameterDFS(G, 0, -1);
        prev.assign((int)G.size(), -1);
        auto [ev, distance2] = DiameterDFS(G, leaf, -1);
        path.clear();
        int cur = ev;
        while (cur != -1) path.push_back(cur), cur = prev[cur];
        return path;
    }
};

// find diameter of weighted tree
template<class Weight, class Graph = vector<vector<pair<int, Weight>>>> struct WeightedDiameter {
    vector<int> path;
    vector<pair<int, Weight>> prev;

    WeightedDiameter() {}
    WeightedDiameter(const Graph &G) {
        solve(G);
    }
    pair<int, Weight> DiameterDFS(const Graph &G, int v, int p) {
        pair<int, Weight> res{v, 0};
        for (auto [to, ew] : G[v]) {
            if (to == p) continue;
            pair<int, Weight> tmp = DiameterDFS(G, to, v);
            tmp.second += ew;
            if (tmp.second > res.second) res = tmp, prev[to] = {v, ew};
        }
        return res;
    }
    pair<Weight, vector<int>> solve(const Graph &G) {
        Weight res = 0;
        prev.assign((int)G.size(), make_pair(-1, -1));
        auto [leaf, distance] = DiameterDFS(G, 0, -1);
        prev.assign((int)G.size(), make_pair(-1, -1));
        auto [ev, distance2] = DiameterDFS(G, leaf, -1);
        path.clear();
        int cur = ev;
        while (cur != -1) {
            if (prev[cur].first != -1) res += prev[cur].second;
            path.push_back(cur), cur = prev[cur].first;
        }
        return {res, path};
    }
};


//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Tree Diameter
void Yosupo_Tree_Diameter() {
    int N, a, b;
    long long c;
    cin >> N;
    vector<vector<pair<int, long long>>> G(N);
    for (int i = 0; i < N-1; i++) {
        cin >> a >> b >> c;
        G[a].emplace_back(b, c), G[b].emplace_back(a, c);
    }
    WeightedDiameter<long long> di;
    auto [len, path] = di.solve(G);
    cout << len << " " << (int)path.size() << endl;
    for (auto v : path) cout << v << " ";
    cout << endl;
}

// AGC 033 C - Removing Coins
void AGC_033_C() {
    int N, a, b;
    cin >> N;
    vector<vector<int>> G(N);
    for (int i = 0; i < N-1; ++i) {
        cin >> a >> b; --a, --b;
        G[a].push_back(b), G[b].push_back(a);
    }
    Diameter di;
    auto path = di.solve(G);    
    if (path.size() % 3 == 2) cout << "Second" << endl;
    else cout << "First" << endl;
}


int main() {
    //Yosupo_Tree_Diameter();
    AGC_033_C();
}