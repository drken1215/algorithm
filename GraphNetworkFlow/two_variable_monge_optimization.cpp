//
// 2 変数 Monge 関数の和の最小化
//
// verified:
//   AtCoder ABC 347 G - Grid Coloring 2 (5 値)
//     https://atcoder.jp/contests/abc347/tasks/abc347_g
//
//   AtCoder ARC 129 E - Yet Another Minimization (M 値)
//     https://atcoder.jp/contests/arc129/tasks/arc129_e
//
//   AtCoder ARC 107 F - Sum of Abs (3 値)
//     https://atcoder.jp/contests/arc107/tasks/arc107_f
//
//   KUPC 2017 H - Make a Potion (バラバラ)
//     https://atcoder.jp/contests/arc107/tasks/arc107_f
//
//   AtCoder ABC 397 G - G - Maximize Distance (段階的な INF 設定)
//     https://atcoder.jp/contests/abc397/tasks/abc397_g
//


#include <bits/stdc++.h>
using namespace std;


// 1, 2, 3-variable submodular optimization 
template<class COST> struct ThreeVariableSubmodularOpt {
    // constructors
    ThreeVariableSubmodularOpt() : N(2), S(0), T(0), OFFSET(0) {}
    ThreeVariableSubmodularOpt(int n, COST inf = numeric_limits<COST>::max() / 2)
    : N(n), S(n), T(n + 1), OFFSET(0), INF(inf), list(n + 2) {}
    
    // initializer
    void init(int n, COST inf = numeric_limits<COST>::max() / 2) {
        N = n, S = n, T = n + 1;
        OFFSET = 0, INF = inf;
        list.clear();
        list.resize(N + 2);
        pos.clear();
    }

    // add constant cost
    void add_cost(COST cost) {
        OFFSET += cost;
    }

    // add 1-variable submodular function
    void add_single_cost(int xi, COST false_cost, COST true_cost) {
        assert(0 <= xi && xi < N);
        if (false_cost >= true_cost) {
            OFFSET += true_cost;
            add_edge(S, xi, false_cost - true_cost);
        } else {
            OFFSET += false_cost;
            add_edge(xi, T, true_cost - false_cost);
        }
    }
    void add_single_cost_01(int xi, COST false_cost, COST true_cost) {
        add_single_cost(xi, false_cost, true_cost);
    }
    
    // add "project selection" constraint
    // xi = T, xj = F: strictly prohibited
    void add_psp_constraint(int xi, int xj) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        add_edge(xi, xj, INF);
    }
    void add_psp_constraint_01(int xi, int xj) {
        add_psp_constraint(xj, xi);
    }
    void add_psp_constraint_10(int xi, int xj) {
        add_psp_constraint(xi, xj);
    }
    
    // add "project selection" penalty
    // xi = T, xj = F: cost C
    void add_psp_penalty(int xi, int xj, COST C) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(C >= 0);
        add_edge(xi, xj, C);
    }
    void add_psp_penalty_01(int xi, int xj, COST C) {
        add_psp_penalty(xj, xi, C);
    }
    void add_psp_penalty_10(int xi, int xj, COST C) {
        add_psp_penalty(xi, xj, C);
    }
    
    // add both True profit
    // xi = T, xj = T: profit P (cost -P)
    void add_both_true_profit(int xi, int xj, COST P) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(P >= 0);
        OFFSET -= P;
        add_edge(S, xi, P);
        add_edge(xi, xj, P);
    }
    
    // add both False profit
    // xi = F, xj = F: profit P (cost -P)
    void add_both_false_profit(int xi, int xj, COST P) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(P >= 0);
        OFFSET -= P;
        add_edge(xj, T, P);
        add_edge(xi, xj, P);
    }
    
    // add general 2-variable submodular function
    // (xi, xj) = (F, F): A, (F, T): B
    // (xi, xj) = (T, F): C, (T, T): D
    void add_submodular_function(int xi, int xj, COST A, COST B, COST C, COST D) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(B + C >= A + D);  // assure submodular function
        OFFSET += A;
        add_single_cost(xi, 0, D - B);
        add_single_cost(xj, 0, B - A);
        add_psp_penalty(xi, xj, B + C - A - D);
    }
    
    // add all True profit
    // y = F: not gain profit (= cost is P), T: gain profit (= cost is 0)
    // y: T, xi: F is prohibited
    void add_all_true_profit(const vector<int> &xs, COST P) {
        assert(P >= 0);
        int y = (int)list.size();
        list.resize(y + 1);
        OFFSET -= P;
        add_edge(S, y, P);
        for (auto xi : xs) {
            assert(xi >= 0 && xi < N);
            add_edge(y, xi, INF);
        }
    }
    
    // add all False profit
    // y = F: gain profit (= cost is 0), T: not gain profit (= cost is P)
    // xi = T, y = F is prohibited
    void add_all_false_profit(const vector<int> &xs, COST P) {
        assert(P >= 0);
        int y = (int)list.size();
        list.resize(y + 1);
        OFFSET -= P;
        add_edge(y, T, P);
        for (auto xi : xs) {
            assert(xi >= 0 && xi < N);
            add_edge(xi, y, INF);
        }
    }
    
    // add general 3-variable submodular function
    // (xi, xj, xk) = (F, F, F): cost A
    // (xi, xj, xk) = (F, F, T): cost B
    // (xi, xj, xk) = (F, T, F): cost C
    // (xi, xj, xk) = (F, T, T): cost D
    // (xi, xj, xk) = (T, F, F): cost E
    // (xi, xj, xk) = (T, F, T): cost F
    // (xi, xj, xk) = (T, T, F): cost G
    // (xi, xj, xk) = (T, T, T): cost H
    void add_submodular_function(int xi, int xj, int xk,
                                 COST A, COST B, COST C, COST D,
                                 COST E, COST F, COST G, COST H) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(0 <= xk && xk < N);
        COST P = (A + D + F + G) - (B + C + E + H);
        COST P12 = (C + E) - (A + G), P13 = (D + G) - (C + H);
        COST P21 = (D + F) - (B + H), P23 = (B + C) - (A + D);
        COST P31 = (B + E) - (A + F), P32 = (F + G) - (E + H);
        assert(P12 >= 0 && P21 >= 0);
        assert(P23 >= 0 && P32 >= 0);
        assert(P31 >= 0 && P13 >= 0);
        if (P >= 0) {
            OFFSET += A;
            add_single_cost(xi, 0, F - B);
            add_single_cost(xj, 0, G - E);
            add_single_cost(xk, 0, D - C);
            add_psp_penalty(xj, xi, P12);
            add_psp_penalty(xk, xj, P23);
            add_psp_penalty(xi, xk, P31);
            add_all_true_profit({xi, xj, xk}, P);
        } else {
            OFFSET += H;
            add_single_cost(xi, C - G, 0);
            add_single_cost(xj, B - D, 0);
            add_single_cost(xk, E - F, 0);
            add_psp_penalty(xi, xj, P21);
            add_psp_penalty(xj, xk, P32);
            add_psp_penalty(xk, xi, P13);
            add_all_false_profit({xi, xj, xk}, -P);
        }
    }
    
    // solve
    COST solve() {
        return dinic() + OFFSET;
    }
    
    // reconstrcut the optimal assignment
    vector<bool> reconstruct() {
        vector<bool> res(N, false), seen(list.size(), false);
        queue<int> que;
        seen[S] = true;
        que.push(S);
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            for (const auto &e : list[v]) {
                if (e.cap && !seen[e.to]) {
                    if (e.to < N) res[e.to] = true;
                    seen[e.to] = true;
                    que.push(e.to);
                }
            }
        }
        return res;
    }
    
    // debug
    friend ostream& operator << (ostream& s, const ThreeVariableSubmodularOpt &tvs) {
        const auto &edges = tvs.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
    
private:
    // edge class
    struct Edge {
        // core members
        int rev, from, to;
        COST cap, icap, flow;
        
        // constructor
        Edge(int r, int f, int t, COST c)
        : rev(r), from(f), to(t), cap(c), icap(c), flow(0) {}
        void reset() { cap = icap, flow = 0; }
        
        // debug
        friend ostream& operator << (ostream& s, const Edge& E) {
            return s << E.from << "->" << E.to << '(' << E.flow << '/' << E.icap << ')';
        }
    };
    
    // inner data
    int N, S, T;
    COST OFFSET, INF;
    vector<vector<Edge>> list;
    vector<pair<int,int>> pos;
    
    // add edge
    Edge &get_rev_edge(const Edge &e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    Edge &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const Edge &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<Edge> get_edges() const {
        vector<Edge> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    void add_edge(int from, int to, COST cap) {
        if (!cap) return;
        pos.emplace_back(from, (int)list[from].size());
        list[from].push_back(Edge((int)list[to].size(), from, to, cap));
        list[to].push_back(Edge((int)list[from].size() - 1, to, from, 0));
    }
    
    // Dinic's algorithm
    COST dinic(COST limit_flow) {
        COST current_flow = 0;
        vector<int> level((int)list.size(), -1), iter((int)list.size(), 0);
        
        // Dinic BFS
        auto bfs = [&]() -> void {
            level.assign((int)list.size(), -1);
            level[S] = 0;
            queue<int> que;
            que.push(S);
            while (!que.empty()) {
                int v = que.front();
                que.pop();
                for (const Edge &e : list[v]) {
                    if (level[e.to] < 0 && e.cap > 0) {
                        level[e.to] = level[v] + 1;
                        if (e.to == T) return;
                        que.push(e.to);
                    }
                }
            }
        };
        
        // Dinic DFS
        auto dfs = [&](auto self, int v, COST up_flow) {
            if (v == T) return up_flow;
            COST res_flow = 0;
            for (int &i = iter[v]; i < (int)list[v].size(); ++i) {
                Edge &e = list[v][i], &re = get_rev_edge(e);
                if (level[v] >= level[e.to] || e.cap == 0) continue;
                COST flow = self(self, e.to, min(up_flow - res_flow, e.cap));
                if (flow <= 0) continue;
                res_flow += flow;
                e.cap -= flow, e.flow += flow;
                re.cap += flow, re.flow -= flow;
                if (res_flow == up_flow) break;
            }
            return res_flow;
        };
        
        // flow
        while (current_flow < limit_flow) {
            bfs();
            if (level[T] < 0) break;
            iter.assign((int)iter.size(), 0);
            while (current_flow < limit_flow) {
                COST flow = dfs(dfs, S, limit_flow - current_flow);
                if (!flow) break;
                current_flow += flow;
            }
        }
        return current_flow;
    };
    COST dinic() {
        return dinic(numeric_limits<COST>::max() / 2);
    }
};

// K-value Two Variable Monge Function Optimization 
/*
    X[i] = 0, 1, ..., K-1 -> (x[i][1], ..., x[i][K-1])

    X[i] < d -> x[i][d] = 1
    X[i] = d -> x[i][1] = 0, ..., x[i][d] = 0, x[i][d+1] = 1, ..., x[i][K] = 1

    X[i] = 0   -> (1, 1, 1, ..., 1, 1)
    X[i] = 1   -> (0, 1, 1, ..., 1, 1)
    X[i] = 2   -> (0, 0, 1, ..., 1, 1)
    ...
    X[i] = K-2 -> (0, 0, 0, ..., 0, 1)
    X[i] = K-1 -> (0, 0, 0, ..., 0, 0)
 */
template<class COST> struct TwoVariableMongeOpt {
    // inner data
    int N, N01;
    COST INF;
    vector<int> ks;  // size of x[i]
    vector<vector<int>> x;  // index of x[i][k] in normal submodular optimization
    ThreeVariableSubmodularOpt<COST> tvs;

    // constructors
    TwoVariableMongeOpt() {}
    TwoVariableMongeOpt(int N, int K, COST inf = numeric_limits<COST>::max() / 2) {
        vector<int> ks(N, K);
        init(ks, inf);
    }
    TwoVariableMongeOpt(const vector<int> &ks, COST inf = numeric_limits<COST>::max() / 2) {
        init(ks, inf);
    }
    void init(const vector<int> &iks, COST inf = numeric_limits<COST>::max() / 2) {
        N = (int)iks.size(), INF = inf, ks = iks, N01 = 0;
        x.resize(N);
        for (int i = 0; i < N; i++) {
            assert(ks[i] >= 2);
            x[i].assign(ks[i], 0);
            for (int k = 1; k < ks[i]; k++) x[i][k] = N01++;
        }
        tvs.init(N01, INF);
        for (int i = 0; i < N; i++) {
            for (int k = 1; k < ks[i] - 1; k++) {
                tvs.add_psp_constraint(x[i][k], x[i][k + 1]);
            }
        }
    }

    // add constant cost
    void add_cost(COST cost) {
        tvs.add_cost(cost);
    }

    // add 1-variable function
    void add_single_cost(int xi, const vector<COST> &cost) {
        assert(0 <= xi && xi < N);
        assert((int)cost.size() == ks[xi]);
        tvs.add_cost(cost[ks[xi] - 1]);
        for (int k = 1; k < ks[xi]; k++) {
            tvs.add_single_cost(x[xi][k], 0, cost[k-1] - cost[k]);
        }
    }

    // add 2-variable Monge function
    // cost[i][j]+cost[i+1][j+1] <= cost[i+1][j]+cost[i][j+1]
    void add_monge_function(int xi, int xj, vector<vector<COST>> cost) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(xi != xj);
        assert(cost.size() == ks[xi]);
        assert(cost[0].size() == ks[xj]);
        vector<COST> icost(ks[xi]), jcost(ks[xj]);
        for (int ki = 0; ki < ks[xi]; ki++) {
            icost[ki] = cost[ki][0];
            for (int kj = 0; kj < ks[xj]; kj++) cost[ki][kj] -= icost[ki];
        }
        for (int kj = 0; kj < ks[xj]; kj++) {
            jcost[kj] = cost[0][kj];
            for (int ki = 0; ki < ks[xi]; ki++) cost[ki][kj] -= jcost[kj];
        }
        add_single_cost(xi, icost), add_single_cost(xj, jcost);
        for (int ki = 1; ki < ks[xi]; ki++) {
            for (int kj = 1; kj < ks[xj]; kj++) {
                COST c = cost[ki][kj] - cost[ki][kj-1] - cost[ki-1][kj] + cost[ki-1][kj-1];
                assert(c <= 0);
                tvs.add_both_false_profit(x[xi][ki], x[xj][kj], -c);
            }
        }
    }

    // solve
    COST solve() {
        return tvs.solve();
    }
    
    // reconstrcut the optimal assignment
    vector<int> reconstruct() {
        vector<int> res(N, 0);
        vector<bool> tres = tvs.reconstruct();
        for (int i = 0; i < N; i++) for (int ki = 1; ki < ks[i]; ki++) {
            res[i] += not tres[x[i][ki]];
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

// AtCoder ABC 347 G - Grid Coloring 2
void ABC_347_G() {
    long long N, INF = 1LL<<45; cin >> N;
    vector<vector<long long>> A(N, vector<long long>(N));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) cin >> A[i][j], A[i][j]--;
    TwoVariableMongeOpt<long long> opt(N * N, 5);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        if (A[i][j] >= 0) {
            vector<long long> cost(5, INF);
            cost[A[i][j]] = 0;
            opt.add_single_cost(i*N+j, cost);
        }
        vector<vector<long long>> cost(5, vector<long long>(5, 0));
        for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) {
            cost[x][y] = (x - y) * (x - y);
        }
        if (i+1 < N) {
            opt.add_monge_function(i*N+j, (i+1)*N+j, cost);
        }
        if (j+1 < N) {
            opt.add_monge_function(i*N+j, i*N+j+1, cost);
        }
    }
    long long res = opt.solve();
    auto x = opt.reconstruct();
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) A[i][j] = x[i*N+j]+1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) cout << A[i][j] << " ";
        cout << endl;
    }
}

// AtCoder ARC 129 E - Yet Another Minimization
#define REP(i, a) for (long long i = 0; i < (long long)(a); i++)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); i++)
void ARC_129_E() {
    long long N, M, INF = 1LL<<45;
    cin >> N >> M;
    vector A(N, vector(M, 0LL)), C(N, vector(M, 0LL)), W(N, vector(N, 0LL));
    REP(i, N) REP(j, M) cin >> A[i][j] >> C[i][j];
    REP(i, N) REP2(j, i+1, N) cin >> W[i][j], W[j][i] = W[i][j];

    TwoVariableMongeOpt<long long> opt(N, M);
    REP(i, N) opt.add_single_cost(i, C[i]);
    REP(i, N) REP2(j, i+1, N) {
        vector cost(M, vector(M, 0LL));
        REP(x, M) REP(y, M) cost[x][y] = W[i][j] * abs(A[i][x] - A[j][y]);
        opt.add_monge_function(i, j, cost);
        
    }
    auto res = opt.solve();
    cout << res << endl;
}

// AtCoder ARC 107 F - Sum of Abs
void ARC_107_F() {
    long long N, M, INF = 1LL<<45;
    cin >> N >> M;
    vector<long long> A(N), B(N), U(M), V(M);
    REP(i, N) cin >> A[i];
    REP(i, N) cin >> B[i];
    REP(i, M) cin >> U[i] >> V[i], U[i]--, V[i]--;

    TwoVariableMongeOpt<long long> opt(N, 3);
    REP(i, N) {
        vector<long long> cost{B[i], A[i], -B[i]};
        opt.add_single_cost(i, cost);
    }
    REP(i, M) {
        vector<vector<long long>> cost = {{0, 0, INF}, {0, 0, 0}, {INF, 0, 0}};
        opt.add_monge_function(U[i], V[i], cost);
    }
    auto res = -opt.solve();
    cout << res << endl;
}

// KUPC 2017 H - Make a Potion
#define ALL(x) x.begin(), x.end()
void KUPC_2017_H() {
    using i128 = __int128_t; 
    long long N, M, INF = 1LL<<60;
    cin >> N >> M;
    vector<long long> V(N), H(N), A(M), X(M), B(M), Y(M);
    REP(i, N) cin >> V[i];
    REP(i, N) cin >> H[i];
    vector<vector<long long>> alts(N);
    REP(i, M) {
        cin >> A[i] >> X[i] >> B[i] >> Y[i], A[i]--, B[i]--;
        alts[A[i]].emplace_back(X[i]);
        if (X[i] > 0) alts[A[i]].emplace_back(X[i]-1);
        alts[B[i]].emplace_back(Y[i]);
        if (Y[i] > 0) alts[B[i]].emplace_back(Y[i]-1);
    }
    vector<int> ks(N);
    REP(i, N) {
        alts[i].emplace_back(0), alts[i].emplace_back(V[i]);
        sort(ALL(alts[i])), alts[i].erase(unique(ALL(alts[i])), alts[i].end());
        ks[i] = alts[i].size();
    }

    TwoVariableMongeOpt<i128> opt(ks);
    REP(i, N) {
        vector<i128> cost(alts[i].size());
        REP(j, alts[i].size()) cost[j] = -H[i] * alts[i][j];
        opt.add_single_cost(i, cost);
    }
    REP(i, M) {
        int a = lower_bound(ALL(alts[A[i]]), X[i]) - alts[A[i]].begin();
        int b = lower_bound(ALL(alts[B[i]]), Y[i]) - alts[B[i]].begin();
        if (A[i] == B[i]) {
            // X[A[i]] >= a かつ X[B[i]] < b を禁止
            vector<i128> cost(alts[A[i]].size(), 0);
            REP2(x, a, b) cost[x] = INF;
            opt.add_single_cost(A[i], cost);
        } else {
            // x[A[i]] >= a かつ X[B[i]] < b を禁止
            vector cost(alts[A[i]].size(), vector<i128>(alts[B[i]].size(), 0));
            REP2(x, a, alts[A[i]].size()) REP(y, b) cost[x][y] = INF;
            opt.add_monge_function(A[i], B[i], cost);
        }
    }
    auto cost = opt.solve();
    cout << -(long long)cost << endl;
}

// AtCoder ABC 397 G - G - Maximize Distance
void ABC_397_G() {
    long long N, M, K, INF = 1LL << 30;
    cin >> N >> M >> K;
    vector<int> U(M), V(M);
    for (int i = 0; i < M; i++) cin >> U[i] >> V[i], U[i]--, V[i]--;

    long long low = -1, high = 40;
    while (high - low > 1) {
        long long d = (high + low) / 2, siz = max(d+1, 2LL);
        TwoVariableMongeOpt<long long> opt(N, siz);

        // 頂点 0 は 0、頂点 N-1 は d
        vector<long long> startcost(siz, INF); startcost[0] = 0;
        vector<long long> goalcost(siz, INF); goalcost[d] = 0;
        opt.add_single_cost(0, startcost);
        opt.add_single_cost(N-1, goalcost);

        // 辺 (U[i], V[i]) について、V[i] - U[i] = 1 でコスト 1、V[i] - U[i] >= 2 は認めない
        vector cost(siz, vector(siz, 0LL));
        for (int i = 0; i <= d; i++) for (int j = i+1; j <= d; j++) cost[i][j] = INF * (j-i-1);
        for (int i = 0; i < d; i++) cost[i][i+1] = 1;
        for (int i = 0; i < M; i++) opt.add_monge_function(U[i], V[i], cost);

        // 解く
        long long optval = opt.solve();
        if (optval <= K) low = d;
        else high = d;
    }
    cout << low << endl;
}


int main() {
    //ABC_347_G();
    //ARC_129_E();
    //ARC_107_F();
    //KUPC_2017_H();
    ABC_397_G();
}