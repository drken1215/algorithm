//
// 2 変数 Monge 関数の和の最小化
//
// verified:
//   AtCoder ABC 326 G - Unlock Achievement (for larger, bigger)
//     https://atcoder.jp/contests/abc326/tasks/abc326_g
//
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


// subroutine: submodular optimization
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
    template<class INT> void add_all_true_profit(const vector<INT> &xs, COST P) {
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
    template<class INT> void add_all_false_profit(const vector<INT> &xs, COST P) {
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

// K-value Two Variable Monge Function Optimization 
/*
    X[i] = 0, 1, ..., K-1 -> (x[i][1], ..., x[i][K-1])
    set X[i] <= d  ⇔  x[i][d] = 1

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
            x[i].assign(ks[i] - 1, 0);
            for (int k = 0; k < ks[i] - 1; k++) x[i][k] = N01++;
        }
        tvs.init(N01, INF);
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < ks[i] - 2; k++) {
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
        for (int k = 0; k < ks[xi] - 1; k++) {
            tvs.add_single_cost(x[xi][k], 0, cost[k] - cost[k + 1]);
        }
    }

    // add 2-variable Monge function
    void add_monge_function(int xi, int xj, const vector<vector<COST>> &cost) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(xi != xj);
        assert((int)cost.size() == ks[xi]);
        assert((int)cost[0].size() == ks[xj]);
        vector<COST> icost(ks[xi], 0), jcost(ks[xj], 0);
        for (int ki = 0; ki < ks[xi]; ki++) icost[ki] = cost[ki][0];
        for (int kj = 1; kj < ks[xj]; kj++) jcost[kj] = cost[ks[xi] - 1][kj] - cost[ks[xi] - 1][0];
        add_single_cost(xi, icost);
        add_single_cost(xj, jcost);
        for (int ki = 0; ki < ks[xi] - 1; ki++) {
            for (int kj = 0; kj < ks[xj] - 1; kj++) {
                COST c = cost[ki][kj + 1] - cost[ki][kj] - cost[ki + 1][kj + 1] + cost[ki + 1][kj];
                assert(c >= 0);
                tvs.add_psp_penalty(x[xi][ki], x[xj][kj], c);
            }
        }
    }

    // add all smaller profit (x[xs[i]] <= a[i])
    template<class INT> void add_all_smaller_profit(const vector<INT> &xs, const vector<INT> &a, COST P) {
        assert(xs.size() == a.size());
        vector<INT> txs;
        for (int i = 0; i < (int)xs.size(); i++) {
            assert(a[i] >= 0);
            if (a[i] >= ks[xs[i]] - 1) continue;
            txs[i].emplace_back(x[xs[i]][a[i]]);  // x <= a equals x[a] = True
        }
        tvs.add_all_true_profit(txs, P);
    }

    // add all larger profit (x[xs[i]] > a[i])
    template<class INT> void add_all_larger_profit(const vector<INT> &xs, const vector<INT> &a, COST P) {
        assert(xs.size() == a.size());
        vector<INT> txs;
        for (int i = 0; i < (int)xs.size(); i++) {
            assert(a[i] < ks[xs[i]] - 1);
            if (a[i] < 0) continue;
            txs.emplace_back(x[xs[i]][a[i]]);  // x > a equals x[a] = False
        }
        tvs.add_all_false_profit(txs, P);
    } 

    // solve
    COST solve() {
        return tvs.solve();
    }
    
    // reconstrcut the optimal assignment
    vector<int> reconstruct() {
        vector<int> res(N, 0);
        vector<bool> tres = tvs.reconstruct();
        for (int i = 0; i < N; i++) for (int ki = 0; ki < ks[i] - 1; ki++) {
            res[i] += not tres[x[i][ki]];
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

#define REP(i, a) for (long long i = 0; i < (long long)(a); i++)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); i++)
#define ALL(x) x.begin(), x.end()

// AtCoder ABC 326 G - Unlock Achievement
// skill level: 0, 1, 2, 3, 4
void ABC_326_G() {
    long long N, M;
    cin >> N >> M;
    vector C(N, 0LL), A(M, 0LL);
    vector X(M, vector(N, 0LL)), L(M, vector(N, 0LL));
    for (int i = 0; i < N; i++) cin >> C[i];
    for (int i = 0; i < M; i++) cin >> A[i];
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            X[i][j] = j;
            cin >> L[i][j], L[i][j]--;
            L[i][j]--;  // x[j] > L[j] でボーナスとなるように
        }
    }
    TwoVariableMongeOpt<long long> opt(N, 5);
    for (int i = 0; i < N; i++) {
        vector<long long> cost(5, 0);
        for (int j = 1; j < 5; j++) cost[j] = cost[j - 1] + C[i];
        opt.add_single_cost(i, cost);
    }
    for (int i = 0; i < M; i++) {
        opt.add_all_larger_profit(X[i], L[i], A[i]);
    }
    long long res = -opt.solve();
    cout << res << endl;
}

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
void ARC_129_E() {
    long long N, M;
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
    ABC_326_G();
    //ABC_347_G();
    //ARC_129_E();
    //ARC_107_F();
    //KUPC_2017_H();
    //ABC_397_G();
}