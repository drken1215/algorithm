//
// maximum stable set, in O(1.381^{n})
//
// cf.
//   wata: 指数時間アルゴリズム入門
//     https://www.slideshare.net/wata_orz/ss-12131479
//
// verified:
//   Yosupo Library Checker - Maximum Independent Set
//     https://judge.yosupo.jp/problem/maximum_independent_set
//
//   CODE THANKS FESTIVAL 2017 G - Mixture Drug
//     https://code-thanks-festival-2017-open.contest.atcoder.jp/tasks/code_thanks_festival_2017_g
//


#include <bits/stdc++.h>
using namespace std;


// find maximum independent set
// graph G must be expressed as adjacent-list
vector<int> maximum_independent_set(const vector<vector<int>> &G) {
    int V = (int)G.size();
    if (V == 1) return {0};
    function<vector<int>(vector<bool>&)> connected_case, general_case;

    connected_case = [&](const vector<bool> &can) -> vector<int> {
        int vma = -1, vmi = -1, ma = -1, mi = V, num = 0;
        vector<int> res;
        vector<vector<int>> G2(V);
        for (int v = 0; v < V; v++) {
            if (!can[v]) continue;
            num++;
            for (auto v2 : G[v]) if (can[v2]) G2[v].emplace_back(v2);
            if (ma < (int)G2[v].size()) ma = (int)G2[v].size(), vma = v;
            if (mi > (int)G2[v].size()) mi = (int)G2[v].size(), vmi = v;
        }
        if (num == 1) return {vmi};
        if (ma <= 2) {
            int v = vmi, p = -1;
            vector<int> path{v};
            vector<bool> ncan(V, false);
            ncan[v] = true;
            for (int i = 0; i < num - 1; i++) {
                int v2 = G2[v][0];
                if (ncan[v2]) v2 = G2[v][1];
                v = v2;
                path.emplace_back(v);
                ncan[v] = true;
            }
            if (mi == 1) {
                for (int i = 0; i < (int)path.size(); i += 2) res.emplace_back(path[i]);
            } else {
                for (int i = 0; i < (int)path.size() - 1; i += 2) res.emplace_back(path[i]);
            }
            return res;
        }
        if (mi < 2) {
            vector<bool> ncan = can;
            ncan[vmi] = false;
            for (auto v : G[vmi]) ncan[v] = false;
            res = general_case(ncan);
            res.emplace_back(vmi);
            return res;
        } else {
            vector<bool> ncan = can;
            ncan[vma] = false;
            res = general_case(ncan);
            ncan = can;
            ncan[vma] = false;
            for (auto v : G[vma]) ncan[v] = false;
            auto tmp = general_case(ncan);
            tmp.emplace_back(vma);
            if (res.size() < tmp.size()) swap(res, tmp);
            return res;
        }
    };

    general_case = [&](vector<bool> &can) -> vector<int> {
        vector<int> res;
        vector<bool> seen(V, false);
        auto dfs = [&](auto &&dfs, int v, vector<bool> &component) -> void {
            seen[v] = true;
            for (auto to : G[v]) if (!seen[to] && can[to]) dfs(dfs, to, component);
            component[v] = true;
        };
        for (int v = 0; v < V; ++v) {
            if (!seen[v] && can[v]) {
                vector<bool> component(V, false);
                dfs(dfs, v, component);
                const auto &tmp = connected_case(component);
                for (auto node : tmp) res.emplace_back(node);
            }
        }
        return res;
    };
    vector<bool> can(V, true);
    return general_case(can);
};


//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Maximum Independent Set
void YosupoMaximumIndependentSet() {
    int N, M;
    cin >> N >> M;
    vector<vector<int>> G(N);
    for (int i = 0; i < M; ++i) {
        int a, b;
        cin >> a >> b;
        G[a].push_back(b), G[b].push_back(a);
    }
    auto res = maximum_independent_set(G);
    cout << res.size() << endl;
    for (auto v : res) cout << v << " ";
    cout << endl;
}


int main () {
    YosupoMaximumIndependentSet();
}