//
// 2 変数劣モジュラ関数のグラフ表現
//
// verified:
//   競プロ典型 90 問 040 - Get More Money（★7）
//     https://atcoder.jp/contests/typical90/tasks/typical90_an
//
//   AtCoder ARC 085 E - MUL (for basid psp)
//     https://atcoder.jp/contests/arc085/tasks/arc085_c
//
//   AtCoder ABC 259 G - Grid Card Game (for basid psp)
//     https://atcoder.jp/contests/abc259/tasks/abc259_g
//
//   AtCoder ARC 176 E - Max Vector (for k-value psp)
//     https://atcoder.jp/contests/arc176/tasks/arc176_e
//
//   AtCoder ABC 326 G - Unlock Achievement (for all-true profit)
//     https://atcoder.jp/contests/abc326/tasks/abc326_g
//
//   AOJ 2903 Board (for general 2-variable submodular function)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2903
//


/*
 N 個の bool 変数 x_0, x_1, ..., x_{N-1} について、以下の形のコストが定められたときの最小コストを求める
 
 ・1 変数 xi に関するコスト (1 変数劣モジュラ関数)
    xi = F のときのコスト, xi = T のときのコスト
 
 ・2 変数 xi, xj 間の関係性についてのコスト (2 変数劣モジュラ関数)
 　　(xi, xj) = (F, F): コスト A
 　　(xi, xj) = (F, T): コスト B
 　　(xi, xj) = (T, F): コスト C
 　　(xi, xj) = (T, T): コスト D
 　(ただし、B + C >= A + D でなければならない)
 
 ・よくある例は、A = B = D = 0, C >= 0 の形である (特に関数化している)
    ・この場合は、特に Project Selection Problem と呼ばれ、ローカルには「燃やす埋める」などとも呼ばれる
    ・xi = T, xj = F のときにコスト C がかかる
 
 ・他に面白い例として、A = B = C = 0, D <= 0 の形もある (これも関数化している)
    ・xi = T, xj = T のときに (-D) の利得が得られる
 */


#include <bits/stdc++.h>
using namespace std;


// 2-variable submodular optimization
template<class COST> struct TwoVariableSubmodularOpt {
    // constructors
    TwoVariableSubmodularOpt() : N(2), S(0), T(0), OFFSET(0) {}
    TwoVariableSubmodularOpt(int n, COST inf = 0)
    : N(n), S(n), T(n + 1), OFFSET(0), INF(inf), list(n + 2) {}
    
    // initializer
    void init(int n, COST inf = 0) {
        N = n, S = n, T = n + 1;
        OFFSET = 0, INF = inf;
        list.assign(N + 2, Edge());
        pos.clear();
    }

    // add 1-Variable submodular functioin
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
    
    // add "project selection" constraint
    // xi = T, xj = F: strictly prohibited
    void add_psp_constraint(int xi, int xj) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        add_edge(xi, xj, INF);
    }
    
    // add "project selection" penalty
    // xi = T, xj = F: cost C
    void add_psp_penalty(int xi, int xj, COST C) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(C >= 0);
        add_edge(xi, xj, C);
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
    friend ostream& operator << (ostream& s, const TwoVariableSubmodularOpt &G) {
        const auto &edges = G.get_edges();
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
        return dinic(numeric_limits<COST>::max());
    }
};



//------------------------------//
// Examples
//------------------------------//

// 競プロ典型 90 問 040 - Get More Money（★7）
void Kyopro_Typical_90_040() {
    // 入力
    int N, W;
    cin >> N >> W;
    vector<int> A(N);
    vector<vector<int>> c(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    for (int i = 0; i < N; ++i) {
        int k;
        cin >> k;
        c[i].resize(k);
        for (int j = 0; j < k; ++j) cin >> c[i][j], --c[i][j];
    }
    
    // 家 i に入らない: F, 家 i に入る: T
    const long long INF = 1LL<<50;
    TwoVariableSubmodularOpt<long long> tvs(N, INF);
    for (int i = 0; i < N; ++i) {
        tvs.add_single_cost(i, 0, W - A[i]);
    }
    
    // 家 v in c[i] に入るためには家 i に入る必要がある
    // つまり、v: T, i: F は禁止
    for (int i = 0; i < N; ++i) {
        for (auto v : c[i]) {
            tvs.add_psp_constraint(v, i);
        }
    }
    cout << -tvs.solve() << endl;
}


// ARC 085 E - MUL
void ARC_085_E() {
    int N;
    cin >> N;
    vector<long long> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    
    // i 個目の宝石を割らない: F, i 個目の宝石を割る: T とする
    const long long INF = 1LL<<55;
    TwoVariableSubmodularOpt<long long> tvs(N, INF);
    for (int i = 0; i < N; ++i) {
        tvs.add_single_cost(i, -a[i], 0);
    }
    
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
            if ((j+1) % (i+1) == 0) {
                // i: T, j: F は禁止
                tvs.add_psp_constraint(i, j);
            }
        }
    }
    cout << -tvs.solve() << endl;
}


// ABC 259 G - Grid Card Game
void ABC_259_G() {
    int H, W;
    cin >> H >> W;
    vector<vector<long long>> A(H, vector<long long>(W));
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) {
        cin >> A[i][j];
        A[i][j] *= -1;
    }
    
    // セットアップ
    const long long INF = 1LL<<50;
    TwoVariableSubmodularOpt<long long> tvs(H + W, INF);
    for (int i = 0; i < H; ++i) {
        long long sum = 0;
        for (int j = 0; j < W; ++j) sum += A[i][j];
        tvs.add_single_cost(i, 0, sum);
    }
    for (int j = 0; j < W; ++j) {
        long long sum = 0;
        for (int i = 0; i < H; ++i) sum += A[i][j];
        tvs.add_single_cost(j+H, sum, 0);
    }
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            if (A[i][j] > 0) tvs.add_psp_constraint(i, j+H);
            else tvs.add_psp_penalty(i, j+H, -A[i][j]);
        }
    }
    cout << -tvs.solve() << endl;
}


// AtCoder ARC 176 E - Max Vector
void ARC_176_E() {
    const int VAL = 500;
    
    int N, M;
    cin >> N >> M;
    vector<int> X(N), Y(N);
    for (int i = 0; i < N; ++i) cin >> X[i];
    for (int i = 0; i < N; ++i) cin >> Y[i];
    vector A(M, vector(N, 0));
    for (int i = 0; i < M; ++i) for (int j = 0; j < N; ++j) cin >> A[i][j];
    
    /*
     j = 0, 1, ..., 499 について
     x[i,j] = 1 (Xi が j より大きい) -> cost 1; 0 (Xi が j 以下である) -> cost 0
     y[i,j] = 0 (Yi が j より大きい) -> cost 1; 1 (Yi が j 以下である) -> cost 0
     z[i] = 1 (X 側を選ぶとき), 0 (Y 側を選ぶとき)
     
     x[i,j] = 0, x[i,j+1] = 1 -> cost INF
     y[i,j] = 1, y[i,j+1] = 0 -> cost INF
     
     x[i,X[i]-1] = 1 である -> 0 に cost INF
     y[i,Y[i]-1] = 0 である -> 1 に cost INF
     
     z[j] = 1 かつ X[i,A[j,i]-1] = 0 -> cost INF
     z[j] = 0 かつ Y[i,A[j,i]-1] = 1 -> cost INF
     */
    const long long INF = 1LL<<50;
    TwoVariableSubmodularOpt<long long> tvs(N * VAL * 2 + M, INF);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < VAL; ++j) {
            tvs.add_single_cost(i * VAL + j, 0, 1);
            tvs.add_single_cost(N * VAL + i * VAL + j, 1, 0);
            if (j + 1 < VAL) {
                tvs.add_psp_constraint(i * VAL + j + 1, i * VAL + j);
                tvs.add_psp_constraint(N * VAL + i * VAL + j, N * VAL + i * VAL + j + 1);
            }
        }
        tvs.add_single_cost(i * VAL + X[i] - 1, INF, 0);
        tvs.add_single_cost(N * VAL + i * VAL + Y[i] - 1, 0, INF);
    }
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < N; ++i) {
            tvs.add_psp_constraint(N * VAL * 2 + j, i * VAL + A[j][i] - 1);
            tvs.add_psp_constraint(N * VAL + i * VAL + A[j][i] - 1, N * VAL * 2 + j);
        }
    }
    cout << tvs.solve() << endl;
}


// ABC 326 G - Unlock Achievement
void ABC_326_G() {
    int N, M;
    cin >> N >> M;
    vector<long long> C(N), A(M);
    vector<vector<long long>> L(M, vector<long long>(N));
    for (int i = 0; i < N; ++i) cin >> C[i];
    for (int i = 0; i < M; ++i) cin >> A[i];
    for (int i = 0; i < M; ++i) for (int j = 0; j < N; ++j) cin >> L[i][j];
    
    // セットアップ
    const long long INF = 1LL<<55;
    TwoVariableSubmodularOpt<long long> tvs(N*4, INF);
    for (int i = 0; i < N*4; ++i) {
        tvs.add_single_cost(i, 0, C[i/4]);
        if (i % 4 != 3) tvs.add_psp_constraint(i+1, i);
    }
    for (int i = 0; i < M; ++i) {
        vector<int> ids;
        for (int j = 0; j < N; ++j) {
            if (L[i][j] > 1) ids.push_back(j*4 + (L[i][j] - 2));
        }
        tvs.add_all_true_profit(ids, A[i]);
    }
    long long res = -tvs.solve();
    cout << res << endl;
}


// AOJ 2093 Board
void AOJ_2903() {
    int n, m;
    cin >> n >> m;
    vector<string> fi(n);
    for (int i = 0; i < n; ++i) cin >> fi[i];
    
    auto get_id = [&](int i, int j) -> int { return i * m + j; };
    
    // 0: 横, 1: 縦
    TwoVariableSubmodularOpt<int> tvs(n * m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (fi[i][j] == '.') continue;
            tvs.add_single_cost(get_id(i, j), 1, 1);
            if (i+1 < n && fi[i+1][j] == '#') {
                // (1, 1) だけ 1 の利得 (-1 のコスト)
                tvs.add_both_true_profit(get_id(i, j), get_id(i+1, j), 1);
            }
            if (j+1 < m && fi[i][j+1] == '#') {
                // (0, 0) だけ 1 の利得 (-1 のコスト)
                tvs.add_both_false_profit(get_id(i, j), get_id(i, j+1), 1);
            }
        }
    }
    cout << tvs.solve() << endl;
}


int main() {
    //Kyopro_Typical_90_040();
    //ARC_085_E();
    //ABC_259_G();
    ARC_176_E();
    //ABC_326_G();
    //AOJ_2903();
}
