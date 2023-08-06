//
// 2-SAT Solver
//
// verified:
//   AtCoder Library Practice Contest H - Two SAT
//     https://atcoder.jp/contests/practice2/tasks/practice2_h
//
//   AtCoder ABC 210 F - Coprime Solitaire
//     https://atcoder.jp/contests/abc210/tasks/abc210_f
//


#include <bits/stdc++.h>
using namespace std;


// Decomposition of Strongly Connected Component
struct SCC {
    using Edge = int;
    using Graph = vector<vector<Edge>>;

    // input
    Graph G, rG;

    // result
    vector<vector<int>> scc;
    vector<int> cmp;
    Graph dag;

    // constructor
    SCC(int N = 0) : G(N), rG(N) {}

    // various methods
    void add_edge(int u, int v) {
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
        dag.assign(scc.size(), vector<Edge>());
        newEdges.clear();
        for (int i = 0; i < (int)G.size(); ++i) {
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
    vector<vector<int>> find_scc(bool to_reconstruct = false) {
        // first dfs
        seen.assign((int)G.size(), false);
        vs.clear();
        for (int v = 0; v < (int)G.size(); ++v) if (!seen[v]) dfs(v);

        // back dfs
        int k = 0;
        scc.clear();
        cmp.assign(G.size(), -1);
        seen.assign(G.size(), false);
        for (int i = (int)G.size()-1; i >= 0; --i) {
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
        add_edge(take_not(x), y);
        add_edge(take_not(y), x);
    }
    
    // main solver (if no solution, return empty)
    vector<int> solve() {
        find_scc();
        for (int i = 0; i < num_variables; ++i) {
            if (cmp[i] == cmp[take_not(i)]) return vector<int>();
            solution[i] = (cmp[i] > cmp[take_not(i)]);
        }
        return solution;
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void ACL_Practice_H() {
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

void ABC_201_F() {
    // 素因数分解
    map<int, set<pair<bool, int>>> primes;
    auto extract_primes = [&](int m, bool which, int id) -> void {
        for (int p = 2; p*p <= m; ++p) if (m % p == 0) {
            primes[p].insert({which, id});
            while (m % p == 0) m /= p;
        }
        if (m > 1) primes[m].insert({which, id});
    };
    
    // 入力 (A: 0, B: 1)
    int N;
    cin >> N;
    vector<int> A(N), B(N);
    for (int i = 0; i < N; ++i) {
        cin >> A[i] >> B[i];
        extract_primes(A[i], false, i);
        extract_primes(B[i], true, i);
    }
    
    // 登録すべきリテラルをメモする
    struct closure {
        bool is1, is2;
        int lit1, lit2;
        closure(bool is1, int lit1, bool is2, int lit2)
        : is1(is1), lit1(lit1), is2(is2), lit2(lit2) {}
    };
    vector<closure> lits;
    int num = N;  // スラックリテラルを含めたリテラルの個数
    for (auto [p, ids] : primes) {
        int num_closures_both = 0;
        for (auto [which, id] : ids) {
            if (which && ids.count({!which, id})) ++num_closures_both;
        }
        
        if (num_closures_both  >= 2) {
            cout << "No" << endl;
            return;
        } else if (num_closures_both == 1) {
            // A[i] も B[i] も p で割り切れるなら、他は p で割れてはいけない
            for (auto [which, id] : ids) {
                if (ids.count({!which, id})) continue;
                lits.push_back(closure(!which, id, !which, id));
            }
        } else {
            int yiter = 0;
            for (auto [which, id] : ids) {
                lits.push_back(closure(!which, id, true, yiter+num));
                if (yiter) {
                    lits.push_back(closure(false, yiter-1+num, true, yiter+num));
                    lits.push_back(closure(false, yiter-1+num, !which, id));
                }
                ++yiter;
            }
            num += yiter;
        }
    }
    
    // 2-SAT を構築
    TwoSATSolver twosat(num);
    for (const auto &cl : lits)
        twosat.add_constraint(cl.is1, cl.lit1, cl.is2, cl.lit2);
    
    // 2-SAT を解く
    const auto &res = twosat.solve();
    cout << (!res.empty() ? "Yes" : "No") << endl;
}


int main() {
    //ACL_Practice_H();
    ABC_201_F();
}
