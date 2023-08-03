//
// Functional Graph の閉路をすべて求める
//
// verified:
//   AtCoder ABC 256 E - Takahashi's Anguish
//     https://atcoder.jp/contests/abc256/tasks/abc256_e
//


#include <bits/stdc++.h>
using namespace std;


// Edge Class
template<class T> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1), val(-1) { }
    Edge(int f, int t, T v = -1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& E) {
        return s << E.from << "->" << E.to;
    }
};

// G[v] := 頂点 v から出ている辺
template<class T> struct CycleDetection {
    // input
    vector<Edge<T>> G;
    
    // intermediate results
    vector<bool> seen, finished;
    vector<int> history;
    
    // constructor
    CycleDetection() { }
    CycleDetection(const vector<Edge<T>> &graph) { init(graph); }
    void init(const vector<Edge<T>> &graph) {
        G = graph;
        seen.assign(G.size(), false);
        finished.assign(G.size(), false);
    }
    
    // return the vertex where cycle is detected
    int search(int v) {
        do {
            seen[v] = true;
            history.push_back(v);
            v = G[v].to;
            if (finished[v]) {
                v = -1;
                break;
            }
        } while (!seen[v]);
        pop_history();
        return v;
    }
    
    // pop history
    void pop_history() {
        while (!history.empty()) {
            int v = history.back();
            finished[v] = true;
            history.pop_back();
        }
    }
    
    // reconstruct
    vector<Edge<T>> reconstruct(int pos) {
        // reconstruct the cycle
        vector<Edge<T>> cycle;
        int v = pos;
        do {
            cycle.push_back(G[v]);
            v = G[v].to;
        } while (v != pos);
        return cycle;
    }
    
    // find cycle, v is the start vertex
    vector<Edge<T>> detect_from_v(int v) {
        int pos = search(v);
        if (pos != -1) return reconstruct(pos);
        else return vector<Edge<T>>();
    }
    
    // find all cycle
    vector<vector<Edge<T>>> detect_all() {
        vector<vector<Edge<T>>> res;
        for (int v = 0; v < (int)G.size(); ++v) {
            if (finished[v]) continue;
            int pos = search(v);
            if (pos == -1) continue;
            const vector<Edge<T>> &cycle = reconstruct(pos);
            if (!cycle.empty()) res.push_back(cycle);
        }
        return res;
    }
};


/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void ABC_256_E() {
    int N;
    cin >> N;
    vector<long long> X(N), C(N);
    for (int i = 0; i < N; ++i) cin >> X[i], --X[i];
    for (int i = 0; i < N; ++i) cin >> C[i];
    
    // Functional Graph 構築
    vector<Edge<long long>> G(N);
    for (int i = 0; i < N; ++i) G[i] = Edge<long long>(i, X[i], C[i]);
    
    // 閉路検出
    using Cycle = vector<Edge<long long>>;
    CycleDetection<long long> cd(G);
    const vector<Cycle> &cycles = cd.detect_all();
    
    // 集計
    long long res = 0;
    for (const auto &cycle : cycles) {
        long long min_cost = 1LL<<60;
        for (const auto &e : cycle) min_cost = min(min_cost, e.val);
        res += min_cost;
    }
    cout << res << endl;
}

int main() {
    ABC_256_E();
}

