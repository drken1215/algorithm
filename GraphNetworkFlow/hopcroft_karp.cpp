//
// Hopcroft-Karp の最大二部マッチング
//
// verified
//   AOJ 1163 カードゲーム
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1163&lang=jp
//

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
using namespace std;


struct HopcroftKarp {
    int sizeL, sizeR;
    vector<vector<int> > list; // left to right
    
    // result
    vector<bool> seen, matched;
    vector<int> level, matching;
    
    // base
    HopcroftKarp(int l, int r) : sizeL(l), sizeR(r), list(l, vector<int>()) { }
    inline vector<int>& operator [] (int i) { return list[i]; }
    inline void addedge(int from, int to) { list[from].push_back(to); }
    inline friend ostream& operator << (ostream& s, const HopcroftKarp& G) {
        s << endl;
        for (int i = 0; i < G.list.size(); ++i) {
            s << i << " : ";
            for (int j = 0; j < G.list[i].size(); ++j) {
                s << G.list[i][j];
                if (j + 1 != G.list[i].size()) s << ", ";
            }
            s << endl;
        }
        return s;
    }
    
    // methods
    void hobfs() {
        queue<int> que;
        for (int left = 0; left < sizeL; ++left) {
            level[left] = -1;
            if (!matched[left]) {
                que.push(left);
                level[left] = 0;
            }
        }
        level[sizeL] = sizeL;
        while (!que.empty()) {
            int left = que.front();
            que.pop();
            for (int i = 0; i < list[left].size(); ++i) {
                int right = list[left][i];
                int next = matching[right];
                if (level[next] == -1) {
                    level[next] = level[left] + 1;
                    que.push(next);
                }
            }
        }
    }
    bool hodfs(int left) {
        if (left == sizeL) return true;
        if (seen[left]) return false;
        seen[left] = true;
        for (int i = 0; i < list[left].size(); ++i) {
            int right = list[left][i];
            int next = matching[right];
            if (level[next] > level[left] && hodfs(next)) {
                matching[right] = left;
                return true;
            }
        }
        return false;
    }
    int solve() {
        seen.assign(sizeL, false);
        matched.assign(sizeL, false);
        level.assign(sizeL+1, -1);
        matching.assign(sizeR, sizeL);
        int res = 0;
        while (true) {
            hobfs();
            seen.assign(sizeL, false);
            bool finished = true;
            for (int left = 0; left < sizeL; ++left) {
                if (!matched[left] && hodfs(left)) {
                    matched[left] = true;
                    ++res;
                    finished = false;
                }
            }
            if (finished) break;
        }
        return res;
    }
};


int GCD(int a, int b) { return b ? GCD(b, a%b) : a; }

int main() {
    int N, M;
    while (cin >> N >> M, N) {
        HopcroftKarp G(N, M);
        vector<int> left(N), right(M);
        for (int i = 0; i < N; ++i) cin >> left[i];
        for (int i = 0; i < M; ++i) cin >> right[i];
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                if (GCD(left[i], right[j]) > 1)
                    G.addedge(i, j);
            }
        }
        cout << G.solve() << endl;
    }
}
