//
// 3 変数劣モジュラ関数のグラフ表現
//
// verified (3 変数は未 verify):
//   競プロ典型 90 問 040 - Get More Money（★7）
//     https://atcoder.jp/contests/typical90/tasks/typical90_an
//
//   AtCoder ARC 085 E - MUL (for basid psp)
//     https://atcoder.jp/contests/arc085/tasks/arc085_c
//
//   AtCoder ABC 259 G - Grid Card Game (for basid psp)
//     https://atcoder.jp/contests/abc259/tasks/abc259_g
//
//   AtCoder ABC 326 G - Unlock Achievement (for all-true profit)
//     https://atcoder.jp/contests/abc326/tasks/abc326_g
//
//   AtCoder ABC 225 G - X (for xi = xj = 1 profit)
//     https://atcoder.jp/contests/abc225/tasks/abc225_g
//
//   AOJ 2903 Board (for general 2-variable submodular function)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2903
//
//   yukicoder No.957 植林
//     https://yukicoder.me/problems/no/957
//


#include <bits/stdc++.h>
using namespace std;


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
    ・この場合は、特に Project Selection Problem と呼ばれ、俗に「燃やす埋める」などとも呼ばれる
    ・xi = T, xj = F のときにコスト C がかかる
 
 ・他に面白い例として、A = B = C = 0, D <= 0 の形もある (これも関数化している)
    ・xi = T, xj = T のときに (-D) の利得が得られる
 
 ・3 変数 xi, xj, xk 間の関係性についてのコスト (3 変数劣モジュラ関数)
 　　(xi, xj, xk) = (F, F, F): コスト A
 　　(xi, xj, xk) = (F, F, T): コスト B
 　　(xi, xj, xk) = (F, T, F): コスト C
 　　(xi, xj, xk) = (F, T, T): コスト D
 　　(xi, xj, xk) = (T, F, F): コスト E
 　　(xi, xj, xk) = (T, F, T): コスト F
 　　(xi, xj, xk) = (T, T, F): コスト G
 　　(xi, xj, xk) = (T, T, T): コスト H
 */

// edge class (for max-flow)
template<class FLOW> struct FlowEdge {
    // core members
    int rev, from, to;
    FLOW cap, icap, flow;
    
    // constructor
    constexpr FlowEdge() noexcept = default;
    constexpr FlowEdge(int rev, int from, int to, FLOW cap, FLOW rcap = 0) 
        : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(rcap) {
    }
    void reset() { 
        flow -= icap - cap;
        cap = icap;
    }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowEdge& e) {
        return s << e.from << " -> " << e.to << " (" << e.cap << ", " << e.flow << ")";
    }
};

// graph class (for max-flow)
template<class FLOW> struct FlowGraph {
    // core members
    vector<vector<FlowEdge<FLOW>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    
    // constructor
    FlowGraph(int n = 0) : list(n) { }
    void init(int n = 0) {
        list.clear(), list.resize(n);
        pos.clear();
    }
    void resize(int n) {
        list.resize(n);
    }
    void clear() {
        list.clear(), pos.clear();
    }
    
    // getter
    vector<FlowEdge<FLOW>> &operator [] (int i) {
        assert(0 <= i && i < (int)list.size());
        return list[i];
    }
    const vector<FlowEdge<FLOW>> &operator [] (int i) const {
        assert(0 <= i && i < (int)list.size());
        return list[i];
    }
    size_t size() const noexcept {
        return list.size();
    }
    size_t size_edegs() const noexcept {
        return pos.size();
    }
    FlowEdge<FLOW> &get_rev_edge(const FlowEdge<FLOW> &e) {
        return list[e.to][e.rev];
    }
    const FlowEdge<FLOW> &get_rev_edge(const FlowEdge<FLOW> &e) const {
        return list[e.to][e.rev];
    }
    FlowEdge<FLOW> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowEdge<FLOW> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowEdge<FLOW>> get_edges() const {
        vector<FlowEdge<FLOW>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() const {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowEdge<FLOW> &e : list[i]) e.reset();
        }
    }
    void change_edge(FlowEdge<FLOW> &e, FLOW new_cap, FLOW new_rcap) {
        assert(new_cap >= 0 && new_rcap >= 0);
        FlowEdge<FLOW> &re = get_rev_edge(e);
        e.cap = new_cap, e.icap = new_cap + new_rcap, e.flow = new_rcap;
        re.cap = new_rcap, re.icap = new_cap + new_rcap, re.flow = new_cap;
    }
    
    // add_edge
    void add_edge(int from, int to, FLOW cap, FLOW rcap = 0) {
        assert(0 <= from && from < (int)list.size() && 0 <= to && to < (int)list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowEdge<FLOW>(to_id, from, to, cap, rcap));
        list[to].push_back(FlowEdge<FLOW>(from_id, to, from, rcap, cap));
    }
    void add_bidirected_edge(int from, int to, FLOW cap) {
        assert(0 <= from && from < (int)list.size() && 0 <= to && to < (int)list.size());
        assert(cap >= 0);
        add_edge(from, to, cap, cap);
    }

    // augment
    FLOW augment(int s, int t, FLOW up_flow = numeric_limits<FLOW>::max()) {
        vector<bool> seen(size(), false);
        auto dfs = [&](auto &&dfs, int v, FLOW up_flow) -> FLOW {
            if (v == t) return up_flow;
            seen[v] = true;
            for (int i = 0; i < (int)list[v].size(); i++) {
                FlowEdge<FLOW> &e = list[v][i], &re = get_rev_edge(e);
                if (seen[e.to] || e.cap <= 0) continue;
                FLOW flow = dfs(dfs, e.to, min(up_flow, e.cap));
                if (flow > 0) {
                    e.cap -= flow, e.flow += flow;
                    re.cap += flow, re.flow -= flow;
                    return flow;
                }
            }  
            return FLOW(0); 
        };
        return dfs(dfs, s, up_flow);
    };

    // find reachable nodes from node s (1: s-domain, -1: t-domain, 0: no reach)
    vector<int> find_cut(int s, int t) const {
        vector<int> res(size(), 0);
        auto dfs_s = [&](auto &&dfs_s, int v) -> void {
            res[v] = 1;
            for (const auto &e : list[v]) {
                if (res[e.to] || e.cap <= 0) continue;
                dfs_s(dfs_s, e.to);
            }
        };
        auto dfs_t = [&](auto &&dfs_t, int v) -> void {
            res[v] = -1;
            for (const auto &e : list[v]) {
                auto re = get_rev_edge(e);
                if (res[e.to] || re.cap <= 0) continue;
                dfs_t(dfs_t, e.to);
            }
        };
        dfs_s(dfs_s, s), dfs_t(dfs_t, t);
        return res;
    }

    // finc cutset
    vector<FlowEdge<FLOW>> find_cutset(int s, int t) const {
        vector<int> cut = find_cut(s, t);
        vector<FlowEdge<FLOW>> res;
        const auto &edges = get_edges();
        for (const auto &e : edges) {
            if (cut[e.from] == 1 && cut[e.to] != 1) {
                res.emplace_back(e);
            }
        }
        return res;
    }

    // check if the s-t flow is feasible
    bool is_feasible(int s, int t) const {
        vector<FLOW> b(list.size(), FLOW(0));
        for (int v = 0; v < (int)list.size(); v++) {
            for (const auto &e : list[v]) {
                b[v] += (e.flow - get_rev_edge(e).flow) / 2;
            }
        }
        if (b[s] + b[t] != 0) return false;
        for (int v = 0; v < (int)list.size(); v++) {
            if (v != s && v != t && b[v] != FLOW(0)) return false;
        }
        return true;
    }
    bool is_feasible(int s, int t, FLOW flow) const {
        vector<FLOW> b(list.size(), FLOW(0));
        for (int v = 0; v < (int)list.size(); v++) {
            for (const auto &e : list[v]) {
                b[v] += (e.flow - get_rev_edge(e).flow) / 2;
            }
        }
        if (b[s] != flow) return false;
        if (b[t] != -flow) return false;
        for (int v = 0; v < (int)list.size(); v++) {
            if (v != s && v != t && b[v] != FLOW(0)) return false;
        }
        return true;
    }

    // decompose flow into s-t simple paths and cycles
    using Path = vector<FlowEdge<FLOW>>;
    pair<vector<Path>, vector<Path>> decompose(int s, int t) const {
        struct Arc {
            int to;
            FLOW rem;
            int eidx;
        };
        assert(is_feasible(s, t));
        vector<vector<Arc>> fg(list.size());
        for (int v = 0; v < (int)list.size(); v++) {
            for (int j = 0; j < (int)list[v].size(); j++) {
                FLOW f = list[v][j].icap - list[v][j].cap;
                if (f > 0) fg[v].push_back({list[v][j].to, f, j});
            }
        }
        vector<int> ptr(list.size(), 0), onpath(list.size(), -1);
        vector<pair<int, int>> route;
        vector<int> used;
        vector<Path> paths, cycles;

        auto next_arc = [&](int v) -> int {
            while (ptr[v] < (int)fg[v].size() && fg[v][ptr[v]].rem <= 0) ptr[v]++;
            return (ptr[v] < (int)fg[v].size() ? ptr[v] : -1);
        };
        auto extract = [&](int begin, bool is_cycle) {
            FLOW mi = numeric_limits<FLOW>::max();
            for (int k = begin; k < (int)route.size(); k++) {
                auto [v, i] = route[k];
                mi = min(mi, fg[v][i].rem);
            }
            vector<FlowEdge<FLOW>> seq;
            for (int k = begin; k < (int)route.size(); k++) {
                auto [v, i] = route[k];
                fg[v][i].rem -= mi;
                FlowEdge<FLOW> e = list[v][fg[v][i].eidx];
                e.flow = mi;
                seq.push_back(e);
            }
            if (is_cycle) cycles.push_back(std::move(seq));
            else paths.push_back(std::move(seq));
        };
        auto walk = [&](int start, bool stop_at_t) {
            route.clear();
            int v = start;
            onpath[v] = 0;
            used.push_back(v);
            while (true) {
                int i = next_arc(v), u = fg[v][i].to;
                route.push_back({v, i});
                if (stop_at_t && u == t) {
                    extract(0, false);
                    break;
                }
                if (onpath[u] != -1) {
                    extract(onpath[u], true);
                    break;
                }
                onpath[u] = (int)route.size();
                used.push_back(u);
                v = u;
            }
            for (int w : used) onpath[w] = -1;
            used.clear();
        };

        // extract all s-t paths
        while (next_arc(s) != -1) walk(s, true);

        // decompose remained circulation into cycles
        for (int v = 0; v < (int)list.size(); v++) while (next_arc(v) != -1) walk(v, false);

        return {paths, cycles};
    }

    // debug
    friend ostream& operator << (ostream& s, const FlowGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// Dinic
template<class FLOW> FLOW Dinic(FlowGraph<FLOW> &G, int s, int t, FLOW limit_flow) {
    assert(0 <= s && s < (int)G.size() && 0 <= t && t < (int)G.size() && s != t);
    FLOW current_flow = 0;
    vector<int> level((int)G.size(), -1), iter((int)G.size(), 0);
    
    // Dinic BFS
    auto bfs = [&]() -> void {
        level.assign((int)G.size(), -1);
        level[s] = 0;
        queue<int> que;
        que.push(s);
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            for (const FlowEdge<FLOW> &e : G[v]) {
                if (level[e.to] < 0 && e.cap > 0) {
                    level[e.to] = level[v] + 1;
                    if (e.to == t) return;
                    que.push(e.to);
                }
            }
        }
    };
    
    // Dinic DFS
    auto dfs = [&](auto self, int v, FLOW up_flow) {
        if (v == t) return up_flow;
        FLOW res_flow = 0;
        for (int &i = iter[v]; i < (int)G[v].size(); ++i) {
            FlowEdge<FLOW> &e = G[v][i], &re = G.get_rev_edge(e);
            if (level[v] >= level[e.to] || e.cap <= 0) continue;
            FLOW flow = self(self, e.to, min(up_flow - res_flow, e.cap));
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
        if (level[t] < 0) break;
        iter.assign((int)iter.size(), 0);
        while (current_flow < limit_flow) {
            FLOW flow = dfs(dfs, s, limit_flow - current_flow);
            if (flow <= 0) break;
            current_flow += flow;
        }
    }
    return current_flow;
};

template<class FLOW> FLOW Dinic(FlowGraph<FLOW> &G, int s, int t) {
    return Dinic(G, s, t, numeric_limits<FLOW>::max());
}

// submodular optimization
template<class COST> struct ThreeVariableSubmodularOpt {
    // Graph
    int N, S, T;
    COST OFFSET, INF;
    FlowGraph<COST> G;

    // constructors
    ThreeVariableSubmodularOpt() : N(2), S(0), T(0), OFFSET(0) {}
    ThreeVariableSubmodularOpt(int n, COST inf = numeric_limits<COST>::max() / 2)
    : N(n), S(n), T(n + 1), OFFSET(0), INF(inf), G(n + 2) {}
    
    // initializer
    void init(int n, COST inf = numeric_limits<COST>::max() / 2) {
        N = n, S = n, T = n + 1;
        OFFSET = 0, INF = inf;
        G.init(N + 2);
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
            if (false_cost - true_cost > 0) G.add_edge(S, xi, false_cost - true_cost);
        } else {
            OFFSET += false_cost;
            G.add_edge(xi, T, true_cost - false_cost);
        }
    }
    void add_single_cost_01(int xi, COST false_cost, COST true_cost) {
        add_single_cost(xi, false_cost, true_cost);
    }
    void add_single_cost_10(int xi, COST false_cost, COST true_cost) {
        add_single_cost(xi, true_cost, false_cost);
    }
    
    // add "project selection" constraint
    // xi = T, xj = F: strictly prohibited
    void add_psp_constraint(int xi, int xj) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(xi != xj);
        G.add_edge(xi, xj, INF);
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
        assert(xi != xj);
        assert(C >= 0);
        if (C > 0) G.add_edge(xi, xj, C);
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
        assert(xi != xj);
        assert(P >= 0);
        OFFSET -= P;
        if (P > 0) G.add_edge(S, xi, P);
        if (P > 0) G.add_edge(xi, xj, P);
    }
    
    // add both False profit
    // xi = F, xj = F: profit P (cost -P)
    void add_both_false_profit(int xi, int xj, COST P) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(xi != xj);
        assert(P >= 0);
        OFFSET -= P;
        if (P > 0) G.add_edge(xj, T, P);
        if (P > 0) G.add_edge(xi, xj, P);
    }
    
    // add general 2-variable submodular function
    // (xi, xj) = (F, F): A, (F, T): B
    // (xi, xj) = (T, F): C, (T, T): D
    void add_submodular_function(int xi, int xj, COST A, COST B, COST C, COST D) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(xi != xj);
        assert(B + C >= A + D);  // assure submodular function
        OFFSET += A;
        add_single_cost(xi, 0, D - B);
        add_single_cost(xj, 0, B - A);
        if (B + C - A - D > 0) add_psp_penalty(xi, xj, B + C - A - D);
    }
    
    // add all True profit
    // y = F: not gain profit (= cost is P), T: gain profit (= cost is 0)
    // y: T, xi: F is prohibited
    void add_all_true_profit(const vector<int> &xs, COST P) {
        assert(P >= 0);
        OFFSET -= P;
        int y = (int)G.size();
        G.resize(y + 1);
        G.add_edge(S, y, P);
        for (auto xi : xs) {
            assert(xi >= 0 && xi < N);
            G.add_edge(y, xi, INF);
        }
    }
    
    // add all False profit
    // y = F: gain profit (= cost is 0), T: not gain profit (= cost is P)
    // xi = T, y = F is prohibited
    void add_all_false_profit(const vector<int> &xs, COST P) {
        assert(P >= 0);
        OFFSET -= P;
        int y = (int)G.size();
        G.resize(y + 1);
        G.add_edge(y, T, P);
        for (auto xi : xs) {
            assert(xi >= 0 && xi < N);
            G.add_edge(xi, y, INF);
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
    COST solve(const string solver = "dinic") {
        if (solver == "dinic") return Dinic(G, S, T) + OFFSET;
        return COST(0);
    }
    
    // reconstrcut the optimal assignment
    vector<bool> reconstruct() {
        vector<bool> res(N, false), seen(G.size(), false);
        queue<int> que;
        seen[S] = true;
        que.push(S);
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            for (const auto &e : G[v]) {
                if (e.cap > 0 && !seen[e.to]) {
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
        const auto &edges = tvs.G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
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
    ThreeVariableSubmodularOpt<long long> tvs(N, INF);
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
    ThreeVariableSubmodularOpt<long long> tvs(N, INF);
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
    ThreeVariableSubmodularOpt<long long> tvs(H + W, INF);
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
    ThreeVariableSubmodularOpt<long long> tvs(N*4, INF);
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


// ABC 225 G - X
void ABC_225_G() {
    long long H, W, C;
    cin >> H >> W >> C;
    vector<vector<long long>> A(H, vector<long long>(W));
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) cin >> A[i][j];
    
    auto get_id = [&](int i, int j) -> int { return i * W + j; };
    
    // セットアップ (F: × を書かない, T: x を書く)
    const long long INF = 1LL<<45;
    ThreeVariableSubmodularOpt<long long> tvs(H * W, INF);
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            tvs.add_single_cost(get_id(i, j), 0, C * 2 - A[i][j]);
            
            // 斜めに隣接すると、C の利得
            if (i+1 < H && j-1 >= 0) {
                tvs.add_both_true_profit(get_id(i, j), get_id(i+1, j-1), C);
            }
            if (i+1 < H && j+1 < W) {
                tvs.add_both_true_profit(get_id(i, j), get_id(i+1, j+1), C);
            }
        }
    }
    
    // 求める
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
    ThreeVariableSubmodularOpt<int> tvs(n * m);
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


// yukicoder No.957 植林
void yukicoder_957() {
    long long H, W;
    cin >> H >> W;
    vector G(H, vector(W, 0LL)); 
    vector R(H, 0LL), C(W, 0LL);
    for (int i = 0; i < H; i++) for (int j = 0; j < W; j++) cin >> G[i][j];
    for (int i = 0; i < H; i++) cin >> R[i];
    for (int j = 0; j < W; j++) cin >> C[j];
    ThreeVariableSubmodularOpt<long long> opt(H + W);
    for (int i = 0; i < H; i++) opt.add_single_cost_10(i, -R[i], 0);
    for (int j = 0; j < W; j++) opt.add_single_cost_10(j+H, -C[j], 0);
    for (int i = 0; i < H; i++) for (int j = 0; j < W; j++) {
        opt.add_submodular_function(i, j+H, 0, G[i][j], G[i][j], G[i][j]);
    }
    cout << -opt.solve() << endl;
}


int main() {
    //Kyopro_Typical_90_040();
    //ARC_085_E();
    //ABC_259_G();
    //ABC_326_G();
    //ABC_225_G();
    //AOJ_2903();
    yukicoder_957();
}