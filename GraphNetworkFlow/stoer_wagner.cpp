//
// 重み付き無向グラフの全域最小カット (by Stoer-Wagner 法, in O(N^3))
//
// verified
//   
//
//


#include <bits/stdc++.h>
using namespace std;


// グラフを隣接行列形式で与える
template<class T> pair<T, vector<int>> stoer_wagner(vector<vector<T>> &G) {
    int V = (int)G.size();
    T best_cost = numeric_limits<T>::max()/2;
    vector<int> best_cut;
    vector<vector<int>> ids(V);
    for (int v = 0; v < V; v++) ids[v] = vector<int>(1, v);
    vector<T> w(V);
    vector<bool> exist(V, true), in_a(V);
    for (int iter = 0; iter < V - 1; iter++) {
        in_a.assign(V, false);
        w.assign(V, 0);
        for (int it = 0, prev; it < V - iter; it++) {
            int sel = -1;
            for (int v = 0; v < V; v++) {
                if (exist[v] && !in_a[v] && (sel == -1 || w[v] > w[sel])) {
                    sel = v;
                }
            }
            if (it == V - iter - 1) {
                if (w[sel] < best_cost) best_cost = w[sel], best_cut = ids[sel];
                ids[prev].insert(ids[prev].end(), ids[sel].begin(), ids[sel].end());
                for (int v = 0; v < V; v++) G[prev][v] = G[v][prev] += G[sel][v];
                exist[sel] = false;
            } else {
                in_a[sel] = true;
                for (int v = 0; v < V; v++) w[v] += G[sel][v];
                prev = sel;
            }
        }
    }
    sort(best_cut.begin(), best_cut.end());
    return {best_cost, best_cut};
}


//------------------------------//
// Examples
//------------------------------//

void user_test() {
    int V;
    cin >> V;
    vector<vector<long long>> G(V, vector<long long>(V, 0));
    for (int i = 0; i < V; i++) for (int j = 0; j < V; j++) cin >> G[i][j];
    auto [cost, cut] = stoer_wagner(G);
    cout << cost << endl;
    for (int i = 0; i < cut.size(); i++) cout << cut[i] << " ";
    cout << endl;
}


int main() {
    user_test();
}