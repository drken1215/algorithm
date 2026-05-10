//
// DAG の最小パス被覆
//
// verified
//   ABC 237 Ex - Hakata
//     https://atcoder.jp/contests/abc237/tasks/abc237_h
//


#include <bits/stdc++.h>
using namespace std;


// DAG min path-cover by Hopcroft-Karp
struct DagPathCover {
    const int NOT_MATCHED = -1;
    
    // input
    int V;
    vector<vector<int>> list; // left to right

    // results
    vector<int> lr, rl;
    
    // intermediate results
    vector<vector<int>> rlist; // right to left
    vector<bool> seen, matched;
    vector<int> level;
    
    // constructor
    DagPathCover(int V) : V(V), list(V), rlist(V) {}
    void add_edge(int from, int to) {
        assert(from >= 0 && from < V);
        assert(to >= 0 && to < V);
        list[from].emplace_back(to);
        rlist[to].emplace_back(from);
    }

    // getter, debugger
    vector<int> &operator [] (int i) { return list[i]; }
    const vector<int> &operator [] (int i) const { return list[i]; }
    constexpr int size() const { return V; }
    friend ostream& operator << (ostream& s, const DagPathCover& G) {
        s << endl;
        for (int i = 0; i < (int)G.list.size(); ++i) {
            s << i << ": ";
            for (int j = 0; j < (int)G.list[i].size(); ++j) {
                s << G.list[i][j];
                if (j + 1 != (int)G.list[i].size()) s << ", ";
            }
            s << endl;
        }
        return s;
    }
    
    // solver
    void hobfs() {
        queue<int> que;
        for (int left = 0; left < V; ++left) {
            level[left] = -1;
            if (!matched[left]) {
                que.push(left);
                level[left] = 0;
            }
        }
        level[V] = V;
        while (!que.empty()) {
            int left = que.front();
            que.pop();
            for (int right : list[left]) {
                int next = rl[right];
                if (next == NOT_MATCHED) next = V;
                if (level[next] == -1) {
                    level[next] = level[left] + 1;
                    que.push(next);
                }
            }
        }
    }
    bool hodfs(int left) {
        if (left == V) return true;
        if (seen[left]) return false;
        seen[left] = true;
        for (int right : list[left]) {
            int next = rl[right];
            if (next == NOT_MATCHED) next = V;
            if (level[next] > level[left] && hodfs(next)) {
                lr[left] = right;
                rl[right] = left;
                return true;
            }
        }
        return false;
    }
    int solve() {
        seen.assign(V, false);
        matched.assign(V, false);
        level.assign(V + 1, -1);
        lr.assign(V, -1);
        rl.assign(V, -1);
        int max_matching = 0;
        while (true) {
            hobfs();
            seen.assign(V, false);
            bool finished = true;
            for (int left = 0; left < V; ++left) {
                if (!matched[left] && hodfs(left)) {
                    matched[left] = true;
                    ++max_matching;
                    finished = false;
                }
            }
            if (finished) break;
        }
        for (int r = 0; r < V; r++) {
            if (rl[r] != NOT_MATCHED) lr[rl[r]] = r;
        }
        return V - max_matching;
    }

    // max stable set
    vector<int> get_stable_set() {
        vector<int> res;
        for (int v = 0; v < V; v++) {
            if (rl[v] == NOT_MATCHED) res.emplace_back(v);
        }
        return res;
    }

    // min path cover
    vector<vector<int>> get_path_cover() {
        auto srcs = get_stable_set();
        vector<vector<int>> res;
        for (auto s : srcs) {
            vector<int> path;
            int v = s;
            while (v != NOT_MATCHED) {
                path.emplace_back(v);
                v = lr[v];
            }
            res.emplace_back(path);
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

// ABC 237 Ex - Hakata
void ABC_237_Ex() {
    string S;
    cin >> S;
    int N = S.size();

    vector<string> ss;
    for (int l = 0; l < N; l++) for (int r = l+1; r <= N; r++) {
        bool ok = true;
        for (int i = l, j = r-1; i < j; i++, j--) if (S[i] != S[j]) ok = false;
        if (ok) ss.emplace_back(S.substr(l, r-l));
    }
    sort(ss.begin(), ss.end());
    ss.erase(unique(ss.begin(), ss.end()), ss.end());
    int V = (int)ss.size();

    // does i include j ?
    auto check = [&](int i, int j) -> bool {
        for (int k = 0; k <= (int)ss[i].size() - (int)ss[j].size(); k++) {
            bool same = true;
            for (int l = 0; l < (int)ss[j].size(); l++) {
                if (ss[i][l+k] != ss[j][l]) same = false;
            }
            if (same) return true;
        }
        return false;
    };

    DagPathCover G(V);
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (i == j) continue;
            if (check(i, j)) G.add_edge(i, j);
        }
    }
    int res = G.solve();
    cout << res << endl;
}


int main() {
    ABC_237_Ex();
}