//
// 二部グラフの最大独立集合, by Hopcroft-Karp, in O(E√V)
//
// verified
//   AtCoder ABC 445 G - Knight Placement
//     https://atcoder.jp/contests/abc445/tasks/abc445_g
//


#include <bits/stdc++.h>
using namespace std;


// Hopcroft-Karp
struct HopcroftKarp {
    const int NOT_MATCHED = -1;
    
    // input
    int size_left, size_right;
    vector<vector<int>> list; // left to right

    // results
    vector<int> lr, rl;
    
    // intermediate results
    vector<bool> seen, matched;
    vector<int> level;
    
    // constructor
    HopcroftKarp(int l, int r) : size_left(l), size_right(r), list(l, vector<int>()) { }
    void add_edge(int from, int to) {
        list[from].push_back(to);
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
};


//------------------------------//
// Examples
//------------------------------//

// AtCoder ABC 445 G - Knight Placement
void ABC_445_G() {
    int N, A, B;
    cin >> N >> A >> B;
    vector dxs{A, A, -A, -A, B, B, -B, -B};
    vector dys{B, -B, B, -B, A, -A, A, -A};
    vector<string> S(N);
    for (int i = 0; i < N; i++) cin >> S[i];

    const int LEFT = 0, RIGHT = 1, NON = -1;
    vector<int> color(N*N, NON);
    auto rec = [&](auto &&rec, int x, int y, int col) -> void {
        color[x*N + y] = col;
        for (int d = 0; d < 8; d++) {
            int x2 = x + dxs[d], y2 = y + dys[d];
            if (x2 < 0 || x2 >= N || y2 < 0 || y2 >= N) continue;
            if (S[x2][y2] == '#') continue;
            if (color[x2*N+y2] == NON) {
                rec(rec, x2, y2, 1-col);
            }
        }
    };
    for (int x = 0; x < N; x++) for (int y = 0; y < N; y++) {
        if (S[x][y] == '#') continue;
        if (color[x*N+y] == NON) rec(rec, x, y, LEFT);
    }

    vector<long long> left, right;
    for (int x = 0; x < N; x++) for (int y = 0; y < N; y++) {
        if (color[x*N+y] == LEFT) left.emplace_back(x*N + y);
        else if (color[x*N+y] == RIGHT) right.emplace_back(x*N + y);
    }
    int L = left.size(), R = right.size();

    HopcroftKarp hk(L, R);
    for (int l = 0; l < L; l++) {
        int v = left[l];
        int x = v / N, y = v % N;
        for (int d = 0; d < 8; d++) {
            int x2 = x + dxs[d], y2 = y + dys[d];
            if (x2 < 0 || x2 >= N || y2 < 0 || y2 >= N) continue;
            if (S[x2][y2] == '#') continue;
            int v2 = x2 * N + y2;
            int r = lower_bound(right.begin(), right.end(), v2) - right.begin();
            hk.add_edge(l, r);
        }
    }
    hk.solve();
    auto res = hk.get_independent_set();

    for (auto [l, id] : res) {
        int v = (l == 0 ? left[id] : right[id]);
        int x = v / N, y = v % N;
        S[x][y] = 'o';
    }
    for (int i = 0; i < N; i++) cout << S[i] << endl;
}


int main() {
    ABC_445_G();
}