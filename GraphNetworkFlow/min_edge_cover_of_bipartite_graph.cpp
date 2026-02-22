//
// 二部グラフの最小辺被覆, by Hopcroft-Karp, in O(E√V)
//
// verified
//
//
//

#include <bits/stdc++.h>
using namespace std;


// Hopcroft-Karp
struct HopcroftKarp {
    const int NOT_MATCHED = -1;
    
    // input
    int size_left, size_right;
    vector<vector<int>> list; // left to right
    vector<vector<int>> rlist; // right to left

    // results
    vector<int> lr, rl;
    
    // intermediate results
    vector<bool> seen, matched;
    vector<int> level;
    
    // constructor
    HopcroftKarp(int L, int R) : size_left(L), size_right(R), list(L), rlist(R) {}
    void add_edge(int from, int to) {
        assert(from >= 0 && from < size_left);
        assert(to >= 0 && to < size_right);
        list[from].emplace_back(to);
        rlist[to].emplace_back(from);
    }

    // getter, debugger
    const vector<int> &operator [] (int i) const { 
        return list[i];
    }
    friend ostream& operator << (ostream& s, const HopcroftKarp& G) {
        s << endl;
        for (int i = 0; i < G.list.size(); ++i) {
            s << i << ": ";
            for (int j = 0; j < G.list[i].size(); ++j) {
                s << G.list[i][j];
                if (j + 1 != G.list[i].size()) s << ", ";
            }
            s << endl;
        }
        return s;
    }
    
    // solver
    void hobfs() {
        queue<int> que;
        for (int left = 0; left < size_left; ++left) {
            level[left] = -1;
            if (!matched[left]) {
                que.push(left);
                level[left] = 0;
            }
        }
        level[size_left] = size_left;
        while (!que.empty()) {
            int left = que.front();
            que.pop();
            for (int i = 0; i < list[left].size(); ++i) {
                int right = list[left][i];
                int next = rl[right];
                if (level[next] == -1) {
                    level[next] = level[left] + 1;
                    que.push(next);
                }
            }
        }
    }
    bool hodfs(int left) {
        if (left == size_left) return true;
        if (seen[left]) return false;
        seen[left] = true;
        for (int i = 0; i < list[left].size(); ++i) {
            int right = list[left][i];
            int next = rl[right];
            if (next == NOT_MATCHED) next = size_left;
            if (level[next] > level[left] && hodfs(next)) {
                rl[right] = left;
                return true;
            }
        }
        return false;
    }
    int solve() {
        seen.assign(size_left, false);
        matched.assign(size_left, false);
        level.assign(size_left + 1, -1);
        lr.assign(size_left, -1);
        rl.assign(size_right, -1);
        int res = 0;
        while (true) {
            hobfs();
            seen.assign(size_left, false);
            bool finished = true;
            for (int left = 0; left < size_left; ++left) {
                if (!matched[left] && hodfs(left)) {
                    matched[left] = true;
                    ++res;
                    finished = false;
                }
            }
            if (finished) break;
        }
        for (int r = 0; r < size_right; r++) {
            if (rl[r] != NOT_MATCHED) lr[rl[r]] = r;
        }
        return res;
    }

    // various construction
    // max matching
    vector<pair<int,int>> get_matching() {
        vector<pair<int,int>> res;
        for (int v = 0; v < size_left; v++) {
            if (lr[v] == NOT_MATCHED) continue;
            res.emplace_back(v, lr[v]);
        }
        return res;
    }

    // enumerate reachable nodes (0: left, 1: right)
    const int LEFT = 0, RIGHT = 1;
    pair<vector<bool>, vector<bool>> get_reachable() {
        vector<bool> can_left(size_left, false);
        vector<bool> can_right(size_right, false);
        queue<pair<int,int>> que;
        for (int v = 0; v < size_left; v++) {
            if (lr[v] == NOT_MATCHED) {
                can_left[v] = true;
                que.push({LEFT, v});
            }
        }
        while (!que.empty()) {
            auto [which, v] = que.front();
            que.pop();
            if (which == LEFT) {
                for (auto r : list[v]) {
                    if (!can_right[r]) {
                        can_right[r] = true;
                        que.push({RIGHT, r});
                    }
                }
            } else {
                int l = rl[v];
                if (l != NOT_MATCHED && !can_left[l]) {
                    can_left[l] = true;
                    que.push({LEFT, l});
                }
            }
        }
        return {can_left, can_right};
    }

    // max independent set (0: left, 1: right)
    vector<pair<int,int>> get_independent_set() {
        vector<pair<int,int>> res;
        auto [can_left, can_right] = get_reachable();
        for (int v = 0; v < size_left; v++) {
            if (can_left[v]) res.emplace_back(LEFT, v);
        }
        for (int v = 0; v < size_right; v++) {
            if (!can_right[v]) res.emplace_back(RIGHT, v);
        }
        return res;
    }

    // min vertex-cover (0: left, 1: right)
    vector<pair<int,int>> get_vertex_cover() {
        vector<pair<int,int>> res;
        auto [can_left, can_right] = get_reachable();
        for (int v = 0; v < size_left; v++) {
            if (!can_left[v]) res.emplace_back(LEFT, v);
        }
        for (int v = 0; v < size_right; v++) {
            if (can_right[v]) res.emplace_back(RIGHT, v);
        }
        return res;
    }

    // min edge-cover (0: left, 1: right)
    vector<pair<int,int>> get_edge_cover() {
        vector<pair<int,int>> res = get_matching();
        for (int v = 0; v < size_left; v++) {
            if (list[v].empty()) return vector<pair<int,int>>();  // infeasible
            if (lr[v] == NOT_MATCHED) res.emplace_back(v, list[v][0]);
        }
        for (int v = 0; v < size_right; v++) {
            if (rlist[v].empty()) return vector<pair<int,int>>();  // infeasible
            if (rl[v] == NOT_MATCHED) res.emplace_back(rlist[v][0], v);
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

int main() {
    
}