//
// 最小費用 b-flow：最小費用循環流に帰着して、Cost-Scaling で解く
//   max-flow 計算を O(f(V, E)) として、O(f(V, E)log(VC))
//
// verified
//   Yosupo Library Checker - Minimum Cost b-flow
//     https://judge.yosupo.jp/problem/min_cost_b_flow
//
//   ABC 421 G - Increase to make it Increasing
//     https://atcoder.jp/contests/abc421/tasks/abc421_g
//
//   KUPC 2014 I - Rain (cost-scaling が苦手説あり: TLE)
//     https://atcoder.jp/contests/kupc2014/tasks/kupc2014_i
//
//   JAG 夏合宿 2013 Day4 I - Multi Path Story (AOJ 2627, cost-scaling が苦手説あり: TLE ギリギリ)
//     https://onlinejudge.u-aizu.ac.jp/challenges/sources/JAG/Summer/2627?year=2013
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


//--------------------------------//
// Max Flow (for feasibility-check)
//--------------------------------//

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


//--------------------------------//
// Min Cost b-Flow
//--------------------------------//

// edge class (for min-cost flow)
template<class FLOW, class COST> struct FlowCostEdge {
    // core members
    int rev, from, to;
    FLOW cap, icap, flow;
    COST cost;
    
    // constructor
    constexpr FlowCostEdge() noexcept = default;
    constexpr FlowCostEdge(int rev, int from, int to, FLOW cap, COST cost)
        : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(0), cost(cost) {
    }
    constexpr FlowCostEdge(int rev, int from, int to, FLOW cap, FLOW rcap, COST cost)
        : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(rcap), cost(cost) {
    }
    void reset() { 
        flow -= icap - cap;
        cap = icap;
    }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowCostEdge& e) {
        return s << e.from << " -> " << e.to << " (" << e.cap << ", " << e.flow << ", " << e.cost << ")";
    }
};

// graph class (for min-cost flow)
template<class FLOW, class COST> struct FlowCostGraph {
    // core members
    vector<vector<FlowCostEdge<FLOW, COST>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    vector<COST> pot;  // pot[v] := potential (e.cost + pot[e.from] - pos[e.to] >= 0)
    bool include_negative_edge = false;
    
    // constructor
    FlowCostGraph(int n = 0) : list(n), pot(n), include_negative_edge(false) { }
    void init(int n = 0) {
        list.clear(), list.resize(n);
        pos.clear();
        pot.assign(n, 0);
        include_negative_edge = false;
    }
    
    // getter
    vector<FlowCostEdge<FLOW, COST>> &operator [] (int i) {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    const vector<FlowCostEdge<FLOW, COST>> &operator [] (int i) const {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowCostEdge<FLOW, COST> &get_rev_edge(const FlowCostEdge<FLOW, COST> &e) {
        return list[e.to][e.rev];
    }
    const FlowCostEdge<FLOW, COST> &get_rev_edge(const FlowCostEdge<FLOW, COST> &e) const {
        return list[e.to][e.rev];
    }
    FlowCostEdge<FLOW, COST> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowCostEdge<FLOW, COST> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowCostEdge<FLOW, COST>> get_edges() const {
        vector<FlowCostEdge<FLOW, COST>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowCostEdge<FLOW, COST> &e : list[i]) e.reset();
        }
    }
    
    // add_edge
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, 0, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, 0, cap, -cost));
        if (cost < 0) include_negative_edge = true;
    }
    void add_edge(int from, int to, FLOW cap, FLOW rcap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, rcap, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, rcap, cap, -cost));
        if (cost < 0) include_negative_edge = true;
    }
    void add_bidirected_edge(int from, int to, FLOW cap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        add_edge(from, to, cap, cap, cost);
    }

    // find initial potential (to resolve initial negative-edge)
    // pot[v] := potential (e.cost + pot[e.from] - pos[e.to] >= 0)
    bool calc_potential_dag() {
        pot.assign(size(), 0);
        vector<int> deg(size(), 0), st;
        for (int v = 0; v < size(); v++) for (const auto &e : list[v]) deg[e.to] += (e.cap > 0);
        st.reserve(size());
        for (int v = 0; v < size(); v++) if (!deg[v]) st.emplace_back(v);
        for (int i = 0; i < size(); i++) {
            if (st.size() == i) return false;  // not DAG
            int cur = st[i];
            for (const auto &e : list[cur]) {
                if (e.cap <= 0) continue;
                deg[e.to]--;
                if (deg[e.to] == 0) st.emplace_back(e.to);
                if (pot[e.to] >= pot[cur] + e.cost) pot[e.to] = pot[cur] + e.cost;
            }
        }
        return true;
    }
    bool calc_potential_spfa() {
        pot.assign(size(), 0);
        queue<int> que;
        vector<bool> inque(size(), false);
        vector<int> cnt(size(), 0);
        for (int v = 0; v < size(); v++) que.push(v), inque[v] = true;
        while (!que.empty()) {
            int cur = que.front();
            que.pop();
            inque[cur] = false;
            if (cnt[cur] > size()) return false;  // include negative-cycle
            cnt[cur]++;
            for (const auto &e : list[cur]) {
                if (e.cap <= 0) continue;
                if (pot[e.to] > pot[cur] + e.cost) {
                    pot[e.to] = pot[cur] + e.cost;
                    if (!inque[e.to]) inque[e.to] = true, que.push(e.to);
                }
            }
        }
        return true;
    }
    bool calc_potential() {
        return calc_potential_dag() || calc_potential_spfa();
    }
    bool init_potential() {
        if (!include_negative_edge) return true;
        return calc_potential();
    }

    // decompose flow into s-t simple paths and cycles
    using Path = vector<FlowCostEdge<FLOW, COST>>;
    pair<vector<Path>, vector<Path>> decompose(int s, int t) const {
        struct Arc {
            int to;
            FLOW rem;
            int eidx;
        };
        vector<vector<Arc>> fg(list.size());
        for (int v = 0; v < (int)list.size(); v++) {
            for (int j = 0; j < (int)list[v].size(); j++) {
                FLOW f = list[v][j].icap - list[v][j].cap;
                if (f > 0) fg[v].push_back({list[v][j].to, f, j});
            }
        }
        vector<Path> paths, cycles;

        auto build = [&](const vector<pair<int,int>> &route, bool is_cycle) {
            FLOW mi = numeric_limits<FLOW>::max();
            for (auto [v,i] : route) mi = min(mi, fg[v][i].rem);
            vector<FlowCostEdge<FLOW,COST>> seq;
            for (auto [v,i] : route) {
                fg[v][i].rem -= mi;
                FlowCostEdge<FLOW,COST> e = list[v][fg[v][i].eidx];
                e.flow = mi;
                seq.push_back(e);
            }
            (is_cycle ? cycles : paths).push_back(move(seq));
        };

        // --- Phase 1: extract all cycles and make graph DAG ---
        vector<int> color(list.size(), 0);          // 0:未訪問 1:スタック上 2:完了
        vector<int> pos_in_stack(list.size(), -1);
        vector<pair<int, int>> stk;
        auto dfs = [&](auto &&dfs, int v) -> bool {
            color[v] = 1; pos_in_stack[v] = (int)stk.size();
            for (int i = 0; i < (int)fg[v].size(); i++) {
                if (fg[v][i].rem <= 0) continue;
                int u = fg[v][i].to;
                if (color[u] == 1) {
                    vector<pair<int,int>> route;
                    for (int k = pos_in_stack[u]; k < (int)stk.size(); k++) {
                        route.push_back(stk[k]);
                    }
                    route.push_back({v, i});
                    build(route, true);
                    return true;
                }
                if (color[u] == 0) {
                    stk.push_back({v, i});
                    if (dfs(dfs, u)) return true;
                    stk.pop_back();
                }
            }
            color[v] = 2; pos_in_stack[v] = -1;
            return false;
        };
        while (true) {
            fill(color.begin(), color.end(), 0);
            stk.clear();
            bool found = false;
            for (int v = 0; v < list.size() && !found; v++) {
                if (color[v] == 0 && dfs(dfs, v)) found = true;
            }
            if (!found) break;
        }

        // --- Phase 2: find all s-t paths ---
        vector<int> ptr(list.size(), 0);
        auto next_arc = [&](int v) -> int {
            while (ptr[v] < (int)fg[v].size() && fg[v][ptr[v]].rem <= 0) ptr[v]++;
            return ptr[v] < (int)fg[v].size() ? ptr[v] : -1;
        };
        while (next_arc(s) != -1) {
            vector<pair<int,int>> route;
            int v = s;
            while (v != t) {
                int i = next_arc(v);
                route.push_back({v, i});
                v = fg[v][i].to;
            }
            build(route, false);
        }
        return {paths, cycles};
    }

    // debug
    friend ostream& operator << (ostream& s, const FlowCostGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// Min Cost Circulation Flow by Cost-Scaling 
template<class FLOW, class COST> COST MinCostCirculation(FlowCostGraph<FLOW, COST> &G) {
    COST eps = 0;
    vector<FLOW> balance(G.size(), 0);
    vector<COST> price(G.size(), 0);
    
    auto newcost = [&](const FlowCostEdge<FLOW, COST> &e) -> COST {
        return e.cost * (COST)G.size() - price[e.from] + price[e.to];
    };

    auto ConstructGaux = [&]() -> void {
        vector<bool> visited(G.size(), false);
        auto dfs = [&](auto &&dfs, int v) -> void {
            visited[v] = true;
            for (int i = 0; i < G[v].size(); ++i) {
                FlowCostEdge<FLOW, COST> &e = G[v][i];
                if (e.cap > 0 && !visited[e.to] && newcost(e) < 0) dfs(dfs, e.to);
            }
        };
        for (int v = 0; v < G.size(); ++v) if (balance[v] > 0) dfs(dfs, v);
        for (int v = 0; v < G.size(); ++v) if (visited[v]) price[v] += eps;
    };

    auto augment_blocking_flow = [&]() -> bool {
        vector<int> iter(G.size(), 0);
        auto augment = [&](auto &&augment, int v, FLOW flow) -> FLOW {
            if (balance[v] < 0) {
                FLOW dif = min(flow, -balance[v]);
                balance[v] += dif;
                return dif;
            }
            for (; iter[v] < G[v].size(); iter[v]++) {
                auto &e = G[v][iter[v]], &re = G.get_rev_edge(e);
                if (e.cap > 0 && newcost(e) < 0) {
                    FLOW dif = augment(augment, e.to, min(flow, e.cap));
                    if (dif > 0) {
                        e.cap -= dif, e.flow += dif;
                        re.cap += dif, re.flow -= dif;
                        return dif;
                    }
                }
            }
            return 0;
        };
        bool finish = true;
        for (int v = 0; v < G.size(); ++v) {
            FLOW flow;
            while (balance[v] > 0 && (flow = augment(augment, v, balance[v])) > 0)
                balance[v] -= flow;
            if (balance[v] > 0) finish = false;
        }
        if (finish) return true;
        else return false;
    };

    for (int v = 0; v < G.size(); ++v) {
        for (int i = 0; i < G[v].size(); ++i) {
            FlowCostEdge<FLOW, COST> &e = G[v][i];
            if (e.cap > 0) eps = max(eps, -e.cost * (COST)G.size());
        }
        price[v] = 0;
    }
    while (eps > 1) {
        eps /= 2;
        for (int v = 0; v < G.size(); ++v) {
            for (int i = 0; i < G[v].size(); ++i) {
                auto &e = G[v][i], &re = G.get_rev_edge(e);
                if (e.cap > 0 && newcost(e) < 0) {
                    FLOW flow = e.cap;
                    balance[e.from] -= flow, balance[e.to] += flow;
                    e.cap -= flow, e.flow += flow;
                    re.cap += flow, re.flow -= flow;
                }
            }
        }
        while (true) {
            ConstructGaux();
            if (augment_blocking_flow()) break;
        }
    }
    COST res = 0;
    const auto &edges = G.get_edges();
    for (const auto &e : edges) res += e.flow * e.cost;
    return res;
}

// Minimum Cost b-flow (come down to min-cost circulation)
template<class FLOW, class COST> struct MinCostBFlow {
    // Edge
    struct Edge {
        int from, to;
        FLOW lower_cap, upper_cap, flow;
        COST cost;

        // debug
        friend ostream& operator << (ostream& s, const Edge& e) {
            return s << e.from << "->" << e.to 
            << '(' << e.lower_cap << '~' << e.upper_cap << ';' << e.cost << ')';
        }
    };

    // inner values
    int V;
    vector<Edge> edges;
    vector<FLOW> lower_dss, upper_dss;  // demand (< 0) and supply (> 0)
    vector<COST> dual;
    FlowCostGraph<FLOW, COST> G;

    // constructor
    MinCostBFlow() {}
    MinCostBFlow(int V) : V(V), lower_dss(V, 0), upper_dss(V, 0) {}

    // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        edges.push_back({from, to, 0, cap, 0, cost});
    }
    void add_edge(int from, int to, FLOW lower_cap, FLOW upper_cap, COST cost) {
        assert(lower_cap <= upper_cap);
        edges.push_back({from, to, lower_cap, upper_cap, 0, cost});
    }
    void set_ds(int v, FLOW ds) {
        assert(0 <= v && v < V);
        lower_dss[v] = ds, upper_dss[v] = ds;
    }
    void set_ds(int v, FLOW lower_ds, FLOW upper_ds) {
        assert(0 <= v && v < V);
        assert(lower_ds <= upper_ds);
        lower_dss[v] = lower_ds, upper_dss[v] = upper_ds;
    }

    // solver
    pair<bool, COST> solve(bool calc_potential = true) {
        // lower_ds, upper_ds -> strict ds
        int super = V;
        vector<FLOW> dss(V + 1, 0);
        for (int i = 0; i < V; i++) {
            if (lower_dss[i] == upper_dss[i]) dss[i] = lower_dss[i];
            else if (lower_dss[i] >= 0) {
                add_edge(super, i, lower_dss[i], upper_dss[i], 0);
            } else if (upper_dss[i] < 0) {
                add_edge(i, super, -upper_dss[i], -lower_dss[i], 0);
            } else {
                add_edge(super, i, upper_dss[i], 0);
                add_edge(i, super, -lower_dss[i], 0);
            }
        }

        // pre-flow lower_cap
        FlowGraph<FLOW> sg(V + 3);
        int s = V + 1, t = V + 2;
        for (const auto &e : edges) {
            dss[e.to] += e.lower_cap, dss[e.from] -= e.lower_cap;
            sg.add_edge(e.from, e.to, e.upper_cap - e.lower_cap);
        }

        // ds -> s, t
        FLOW ssum = 0, tsum = 0;
        for (int i = 0; i < V + 1; i++) {
            if (dss[i] > 0) ssum += dss[i], sg.add_edge(s, i, dss[i]);
            else if (dss[i] < 0) tsum -= dss[i], sg.add_edge(i, t, -dss[i]);
        }

        // feasibility check
        if (ssum != tsum) return {false, COST(0)};
        if (Dinic(sg, s, t) < ssum) return {false, COST(0)};

        // come down to min-cost circulation
        G.init(V + 1);
        for (int i = 0; i < (int)edges.size(); i++) {
            auto &e = edges[i];
            const auto &ge = sg.get_edge(i);
            G.add_edge(ge.from, ge.to, ge.cap, ge.flow, e.cost);
        }
        MinCostCirculation(G);

        // find min-cost
        COST res = 0;
        for (int i = 0; i < (int)edges.size(); i++) {
            auto &e = edges[i];
            const auto &ge = G.get_edge(i);
            e.flow = e.upper_cap - ge.cap;
            res += e.flow * e.cost;
        }

        // find dual
        if (calc_potential) {
            G.calc_potential();
            dual = G.pot;
            dual.pop_back();  // eliminate super-node
        }
        return {true, res};
    }
};


//------------------------------//
// Solver
//------------------------------//

// Yosupo Libray Checker - Minimum Cost b-flow
// int 128
using i128 = __int128_t;
i128 to_integer(const string &s) {
    i128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    if (s[0] == '-') res *= -1;
    return res;
}
istream& operator >> (istream &is, i128 &x) {
    string s;
    is >> s;
    x = to_integer(s);
    return is;
}
ostream& operator << (ostream &os, const i128 &x) {
    i128 ax = (x >= 0 ? x : -x);
    char buffer[128];
    char *d = end(buffer);
    do {
         --d;
        *d = "0123456789"[ax % 10];
        ax /= 10;
    } while (ax != 0);
    if (x < 0) {
        --d;
        *d = '-';
    }
    int len = end(buffer) - d;
    if (os.rdbuf()->sputn(d, len) != len) {
        os.setstate(ios_base::badbit);
    }
    return os;
}
void Yosupo_Minimum_Cost_b_flow() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);

    int N, M;
    cin >> N >> M;
    MinCostBFlow<long long, i128> mcf(N);
    vector<i128> B(N);
    for (int i = 0; i < N; i++) cin >> B[i], mcf.set_ds(i, B[i]);
    vector<int> s(M), t(M);
    vector<i128> l(M), u(M), c(M);
    for (int i = 0; i < M; i++) {
        cin >> s[i] >> t[i] >> l[i] >> u[i] >> c[i];
        mcf.add_edge(s[i], t[i], l[i], u[i], c[i]);
    }

    auto [exist, res] = mcf.solve();
    if (!exist) cout << "infeasible" << '\n';
    else {
        cout << res << '\n';
        for (int i = 0; i < N; i++) cout << mcf.dual[i] << '\n';
        for (int i = 0; i < M; i++) cout << mcf.edges[i].flow << '\n';
    }
}

// ABC 421 G - Increase to make it Increasing
void ABC_421_G() {
    long long N, M, INF = 1LL<<45; cin >> N >> M;
    vector<long long> A(N), D(N+1, INF);
    D[0] = 0;
    for (int i = 0; i < N; i++) {
        cin >> A[i];
        if (i) D[i] = A[i] - A[i-1];
    }

    MinCostBFlow<long long, long long> G(N+1);
    for (int v = 0; v <= N; v++) {
        if (D[v] >= 0) G.set_ds(v, 0, D[v]);
        else G.set_ds(v, -INF, D[v]);
    }
    for (int i = 0; i < M; i++) {
        long long u, v; cin >> u >> v; u--;
        G.add_edge(v, u, INF, 1);
    }
    auto [flag, cost] = G.solve();
    cout << (flag ? cost : -1) << endl;
}

// KUPC 2014 I - Rain
void KUPC_2014_I() {
    int N, M, K, INF = 1LL<<20;
    cin >> N >> M >> K;
    vector<long long> C(K), A(M), B(M), D(M), b(N, 0); 
    for (int i = 0; i < K; i++) cin >> C[i], C[i]--;
    for (int i = 0; i < M; i++) cin >> A[i] >> B[i] >> D[i], A[i]--, B[i]--;
    for (int i = 0; i < K; i++) b[A[C[i]]]++, b[B[C[i]]]--;
    MinCostBFlow<int, long long> G(N);
    for (int i = 0; i < N; i++) G.set_ds(i, b[i]);
    for (int i = 0; i < M; i++) G.add_edge(B[i], A[i], INF, D[i]);
    auto [flag, cost] = G.solve();
    cout << (flag ? cost : -1) << endl;
}

// JAG 夏合宿 2013 Day4 I - Multi Path Story (AOJ 2627)
void AOJ_2627() {
    long long N, INF = 1LL << 45;
    cin >> N;
    MinCostBFlow<long long, long long> G(N + 1);
    long long t = N;
    for (int v = 0; v < N; v++) {
        int D;
        cin >> D;
        for (int i = 0; i < D; i++) {
            long long to, w;
            cin >> to >> w, to--;
            G.add_edge(v, to, 1, INF, w);
        }
    }
    for (int v = 1; v < N; v++) G.add_edge(v, t, 0, INF, 0);
    G.add_edge(t, 0, 0, INF, 0);
    auto [flag, cost] = G.solve();
    cout << cost << endl;
}


int main() {
    //Yosupo_Minimum_Cost_b_flow();
    //ABC_421_G();
    //KUPC_2014_I();
    AOJ_2627();
}