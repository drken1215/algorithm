//
// ハンガリアン法 (L <= R を仮定, in O(L^2R))
//    双対変数（次の双対問題の最適解）も求める
//       ・L = R のとき
//           maximize: sum(x_i) + sum(y_j) subject to x_i + y_j <= C_{i,j}
//       ・L < R のとき
//           maximize: sum(x_i) + sum(y_j) subject to x_i + y_j <= C_{i,j}, y_j <= 0
//
// Example
//   Library Checker - Assignment Problem
//     https://judge.yosupo.jp/problem/assignment
//


#include <bits/stdc++.h>
using namespace std;


// Hungarian
template<class VAL> struct Hungarian {
    // inner values
    int L, R;  // left size, right size
    vector<vector<VAL>> G;

    // results
    vector<int> lr, rl;
    vector<VAL> dualL, dualR;
    
    // constructor
    explicit Hungarian(const vector<vector<VAL>> &G_) 
        : L((int)G_.size()), R((int)G_[0].size()), G(G_) { 
        assert(L <= R);
    }

    // getter
    vector<int> get_lr() const { return lr; };
    vector<int> get_rl() const { return rl; };
    pair<vector<VAL>, vector<VAL>> get_dual() const { return {dualL, dualR}; };
    
    // solver
    VAL solve() {
        lr.assign(L, -1), rl.assign(R, -1);
        dualL.assign(L, VAL(0)), dualR.assign(R, VAL(0));
        vector<VAL> dist(R);
        vector<int> index(R), prev(R);
        iota(index.begin(), index.end(), 0);

        auto calc_residue = [&](int i, int j) -> VAL { return G[i][j] - dualR[j]; };
        for (int f = 0; f < L; f++) {
            for (int j = 0; j < R; j++) dist[j] = calc_residue(f, j), prev[j] = f;
            VAL w = 0;
            int j = 0, l = 0;
            for (int s = 0, t = 0;;) {
                if (s == t) {
                    l = s, w = dist[index[t++]];
                    for (int k = t; k < R; k++) {
                        j = index[k];
                        if (dist[j] <= w) {
                            if (dist[j] < w) t = s, w = dist[j];
                            index[k] = index[t], index[t++] = j;
                        }
                    }
                    for (int k = s; k < t; k++) {
                        j = index[k];
                        if (rl[j] < 0) goto augment;
                    }
                }
                int q = index[s++], i = rl[q];
                for (int k = t; k < R; k++) {
                    j = index[k];
                    VAL h = calc_residue(i, j) - calc_residue(i, q) + w;
                    if (h < dist[j]) {
                        dist[j] = h, prev[j] = i;
                        if (h == w) {
                            if (rl[j] < 0) goto augment;
                            index[k] = index[t], index[t++] = j;
                        }
                    }
                }
            }
        augment:
            for (int k = 0; k < l; k++) dualR[index[k]] += dist[index[k]] - w;
            int i = 0;
            do {
                rl[j] = i = prev[j];
                swap(j, lr[i]);
            } while (i != f);
        }
        VAL res = 0;
        for (int i = 0; i < L; i++) {
            res += G[i][lr[i]];
            dualL[i] = G[i][lr[i]] - dualR[lr[i]];
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

// Library Checker - Assignment Problem
void LibraryChecker_Hungarian() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N;
    cin >> N;
    vector<vector<long long>> A(N, vector<long long>(N));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) cin >> A[i][j];
    Hungarian<long long> hung(A);
    long long res = hung.solve();
    auto p = hung.get_lr();
    cout << res << endl;
    for (int i = 0; i < (int)p.size(); i++) cout << p[i] << " ";
    cout << endl;
}


int main() {
    LibraryChecker_Hungarian(); 
}