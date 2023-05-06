//
// 2-SAT Solver
//
// verified:
//   AtCoder Library Practice Contest H - Two SAT
//     https://atcoder.jp/contests/practice2/tasks/practice2_h
//


#include <bits/stdc++.h>
using namespace std;


// Decomposition of Strongly Connected Component
struct SCC {
    using Edge = int;
    using SGraph = vector<vector<Edge>>;

    // input
    SGraph G, rG;

    // result
    vector<vector<int>> scc;
    vector<int> cmp;
    SGraph dag;

    // constructor
    SCC(int N = 0) : G(N), rG(N) {}

    // add edge
    void addedge(int u, int v) {
        G[u].push_back(v);
        rG[v].push_back(u);
    }

    // decomp
    vector<bool> seen;
    vector<int> vs, rvs;
    void dfs(int v) {
        seen[v] = true;
        for (auto e : G[v]) if (!seen[e]) dfs(e);
        vs.push_back(v);
    }
    void rdfs(int v, int k) {
        seen[v] = true;
        cmp[v] = k;
        for (auto e : rG[v]) if (!seen[e]) rdfs(e, k);
        rvs.push_back(v);
    }

    // reconstruct
    set<pair<int,int>> newEdges;
    void reconstruct() {
        int N = (int)G.size();
        int dV = (int)scc.size();
        dag.assign(dV, vector<Edge>());
        newEdges.clear();
        for (int i = 0; i < N; ++i) {
            int u = cmp[i];
            for (auto e : G[i]) {
                int v = cmp[e];
                if (u == v) continue;
                if (!newEdges.count({u, v})) {
                    dag[u].push_back(v);
                    newEdges.insert({u, v});
                }
            }
        }
    }

    // main solver
    vector<vector<int>> find_scc(bool to_reconstruct = true) {
        // first dfs
        int N = (int)G.size();
        seen.assign(N, false);
        vs.clear();
        for (int v = 0; v < N; ++v) if (!seen[v]) dfs(v);

        // back dfs
        int k = 0;
        scc.clear();
        cmp.assign(N, -1);
        seen.assign(N, false);
        for (int i = N - 1; i >= 0; --i) {
            if (!seen[vs[i]]) {
                rvs.clear();
                rdfs(vs[i], k++);
                scc.push_back(rvs);
            }
        }

        // reconstruct DAG
        if (to_reconstruct) reconstruct();
        return scc;
        
    }
};

// 2-SAT Solver
struct TwoSATSolver : SCC {
    // input
    int num_variables;
    
    // result
    vector<int> solution;
    
    // constructor
    TwoSATSolver(int N = 0) : SCC(N * 2), num_variables(N), solution(N) {}
    
    // not
    inline int take_not(int x) {
        if (x < num_variables) return x + num_variables;
        else return x - num_variables;
    }
    
    // add closure
    void add_constraint(bool is_x_true, int x, bool is_y_true, int y) {
        assert(x >= 0 && x < num_variables);
        assert(y >= 0 && y < num_variables);
        if (!is_x_true) x = take_not(x);
        if (!is_y_true) y = take_not(y);
        addedge(take_not(x), y);
        addedge(take_not(y), x);
    }
    
    // main solver
    vector<int> solve() {
        find_scc();
        for (int i = 0; i < num_variables; ++i) {
            // no solution
            if (cmp[i] == cmp[take_not(i)]) {
                return vector<int>();
            }
            solution[i] = (cmp[i] > cmp[take_not(i)]);
        }
        return solution;
    }
};


int main() {
    int N, D;
    cin >> N >> D;
    vector<int> X(N), Y(N);
    for (int i = 0; i < N; ++i) cin >> X[i] >> Y[i];
    
    // 2-SAT Solver
    TwoSATSolver twosat(N);
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
            if (abs(X[i] - X[j]) < D)
                twosat.add_constraint(false, i, false, j);
            if (abs(X[i] - Y[j]) < D)
                twosat.add_constraint(false, i, true, j);
            if (abs(Y[i] - X[j]) < D)
                twosat.add_constraint(true, i, false, j);
            if (abs(Y[i] - Y[j]) < D)
                twosat.add_constraint(true, i, true, j);
        }
    }
    const auto &res = twosat.solve();
                
    if (res.empty()) cout << "No" << endl;
    else {
        cout << "Yes" << endl;
        for (int i = 0; i < N; ++i) {
            cout << (res[i] ? X[i] : Y[i]) << endl;
        }
    }
}
