#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/problems/2171"
#define ERROR 1e-8

#include <bits/stdc++.h>
using namespace std;

#include "double_matrix.hpp"

int main() {
    int N, s, t;
    while (cin >> N >> s >> t, N) {
        --s, --t;
        vector<int> q(N);
        vector<vector<int>> a(N, vector<int>(N));
        for (int i = 0; i < N; ++i) cin >> q[i];
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) cin >> a[i][j];

        // Dijkstra
        const int INF = 1<<29;
        vector<int> dist(N, INF);
        vector<bool> seen(N, 0);
        dist[t] = 0;
        for (int iter = 0; iter < N; ++iter) {
            int curd = INF;
            int v = -1;
            for (int i = 0; i < N; ++i) {
                if (seen[i]) continue;
                if (curd > dist[i]) {
                    curd = dist[i];
                    v = i;
                }
            }
            if (v == -1) break;
            for (int w = 0; w < N; ++w) {
                if (w == v) continue;
                if (a[v][w] == 0) continue;
                dist[w] = min(dist[w], curd + a[v][w]);
            }
            seen[v] = true;
        }
        if (dist[s] >= INF) {
            cout << "impossible" << endl;
            continue;
        }

        // 連立一次方程式を作る
        DoubleMatrix<double> A(N, N, 0);
        vector<double> b(N, 0);
        for (int v = 0; v < N; ++v) {
            if (v == t) {
                A[v][v] = 1;
                b[v] = 0;
            }
            else {
                vector<int> neigbor;
                for (int w = 0; w < N; ++w) {
                    if (a[v][w] == 0) continue;
                    if (q[v] == 1 && dist[w] + a[v][w] != dist[v]) continue;
                    neigbor.push_back(w);
                }
                int K = neigbor.size();
                for (auto w : neigbor) {
                    A[v][w] -= 1;
                    b[v] += a[v][w];
                }
                A[v][v] += K;
            }
        }
          
        // 解く
        vector<double> res;
        auto rank = linear_equation(A, b, res);
        if (res.empty()) cout << "impossible" << endl;
        else cout << fixed << setprecision(15) << res[s] << endl;
    }
}
