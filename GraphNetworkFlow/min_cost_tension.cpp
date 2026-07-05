//
// 最小費用テンション問題
//
// reference:
//   双対性
//     https://www.slideshare.net/slideshow/ss-91375739/91375739
//
// verified
//   JAG 夏合宿 2010 Day4 I - How to Create a Good Game (AOJ 2230)
//     https://onlinejudge.u-aizu.ac.jp/problems/2230
//
//   AtCoder ABC 347 G - Grid Coloring 2
//     https://atcoder.jp/contests/abc347/tasks/abc347_g
//
//   AtCoder ABC 397 G - Maximize Distance
//     https://atcoder.jp/contests/abc397/tasks/abc397_g
//
//   HUPC 2020 H - Traditional Company (AOJ 3171)
//     https://onlinejudge.u-aizu.ac.jp/problems/3171
//
//   AtCoder ABC 393 G - Unevenness (for using frac<i128>)
//     https://atcoder.jp/contests/abc393/tasks/abc393_g
//

/*
    min_{p}: 
        Σ_{v} b(v)p(v) + Σ_{e} {c(e) max(0, p(v) - p(u) - l(e)}
    s.t.
        p(v) - p(u) <= d(e)
    ->
        b-flow (with demand-suply: b)
        edge 
            obj func: e = (u, v) with capacity c(e), cost l(e)
            constraint: e = (u, v) with capacity INF, cost d(e)
        optimal value *= -1

    in general:
        Σ_{v} b(v)p(v) + Σ_{e} f(p(v) - p(u)), where f is concave
*/


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


//--------------------------------//
// Min-cost b-flow (subroutine)
//--------------------------------//

// edge class (for max-flow)
template<class FLOW> struct FlowEdge {
    // core members
    int rev, from, to;
    FLOW cap, icap, flow;
    
    // constructor
    FlowEdge() {}
    FlowEdge(int rev, int from, int to, FLOW cap, FLOW rcap = 0) 
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
        assert(0 <= i && i < list.size());
        return list[i];
    }
    const vector<FlowEdge<FLOW>> &operator [] (int i) const {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowEdge<FLOW> &get_rev_edge(const FlowEdge<FLOW> &e) {
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
    void reset() {
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
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowEdge<FLOW>(to_id, from, to, cap, rcap));
        list[to].push_back(FlowEdge<FLOW>(from_id, to, from, rcap, cap));
    }
    void add_bidirected_edge(int from, int to, FLOW cap) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
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
    vector<int> find_cut(int s, int t) {
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
    bool is_feasible(int s, int t) {
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
    bool is_feasible(int s, int t, FLOW flow) {
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

    // debug
    friend ostream& operator << (ostream& s, const FlowGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// Dinic
template<class FLOW> FLOW Dinic(FlowGraph<FLOW> &G, int s, int t, FLOW limit_flow) {
    assert(0 <= s && s < G.size() && 0 <= t && t < G.size() && s != t);
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

// edge class (for min-cost flow)
template<class FLOW, class COST> struct FlowCostEdge {
    // core members
    int rev, from, to;
    FLOW cap, icap, flow;
    COST cost;
    
    // constructor
    FlowCostEdge() {}
    FlowCostEdge(int rev, int from, int to, FLOW cap, COST cost)
        : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(0), cost(cost) {
    }
    FlowCostEdge(int rev, int from, int to, FLOW cap, FLOW rcap, COST cost)
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


//--------------------------------//
// Min-cost tension
//--------------------------------//

// min-cost tension
/*
    min_{p}: 
        Σ_{v} b(v)p(v) + Σ_{e} {c(e) max(0, p(v) - p(u) - l(e)}
    s.t.
        p(v) - p(u) <= d(e)
    ->
        b-flow (with demand-suply: b)
        edge 
            obj func: e = (u, v) with capacity c(e), cost l(e)
            constraint: e = (u, v) with capacity INF, cost d(e)
        optimal value *= -1

    in general:
        Σ_{v} b(v)p(v) + Σ_{e} f(p(v) - p(u)), where f is concave
*/
template<class FLOW, class COST> struct MinCostTension {
    // inner values
    int N;
    COST OFFSET = 0;
    MinCostBFlow<FLOW, COST> opt;

    // constructor
    MinCostTension() : OFFSET(0) {}
    MinCostTension(int n) : N(n), OFFSET(0), opt(N) {}
    void init(int n) {
        N = n;
        OFFSET = 0;
        opt.init(n);
    }

    // add constant cost
    void add_cost(COST cost) {
        OFFSET += cost;
    }

    // add the part of obj func Σ_{v}b(v)p(v)
    void add_single_coef(int v, FLOW b) {
        assert(0 <= v && v < N);
        assert(opt.lower_dss[v] == opt.upper_dss[v]);
        opt.set_ds(v, opt.lower_dss[v] + b);
    }

    // add tha part of obj func Σ_{e} {c(e) max(0, p(v) - p(u) - l(e)}
    void add_tension_cost(int u, int v, FLOW c, COST l) {
        assert(0 <= u && u < N);
        assert(0 <= v && v < N);
        assert(u != v);
        assert(c >= 0);
        opt.add_edge(u, v, c, l);
    }

    // add constraint p(v) - p(u) <= d
    void add_tension_constraint(int u, int v, FLOW inf, COST d) {
        assert(0 <= u && u < N);
        assert(0 <= v && v < N);
        assert(u != v);
        opt.add_edge(u, v, inf, d);
    }

    // テンション p[v] - p[u] に関する区分線形凸関数 f を足す
    // f を (min_{f}, 傾きが 0 以下・0 以上の部分の傾きの変化点の多重集合）で表す
    // 変化点の多重集合を (変化点, 変化量) の vector で表す
    void add_tension_convex_function(int u, int v, 
    COST mif, const vector<pair<COST,FLOW>> &left, const vector<pair<COST,FLOW>> &right) {
        assert(0 <= u && u < N);
        assert(0 <= v && v < N);
        assert(u != v);
        add_cost(mif);
        for (auto [x, d] : left) add_tension_cost(v, u, d, -x);
        for (auto [x, d] : right) add_tension_cost(u, v, d, x);
    }

    // solver
    pair<bool, COST> solve(bool calc_potential = true) {
        auto [flag, cost] = opt.solve(calc_potential);
        return make_pair(flag, OFFSET - cost);
    }
    vector<FLOW> reconstruct() {
        return opt.dual;
    }
};


//------------------------------//
// Solver
//------------------------------//

// JAG 夏合宿 2010 Day4 I - How to Create a Good Game (AOJ 2230)
/*
    最長路の長さを D として
    max: Σ_{e}(p(v) - p(u) - w(e)) -> min: Σ_{e}w(e) + Σ_{v}outdeg(v)p(v) - Σ_{v}indeg(v)p(v)
    s.t. 
        p(N-1) - p(0) <= D
        p(v) - p(u) >= w(e)
*/
void AOJ_2230() {
    int N, M, sum = 0;
    cin >> N >> M;
    vector<vector<pair<int,int>>> dag(N);
    vector<int> U(M), V(M), W(M), outdeg(N, 0), indeg(N, 0);
    for (int i = 0; i < M; i++) {
        cin >> U[i] >> V[i] >> W[i];
        sum += W[i];
        dag[U[i]].emplace_back(V[i], W[i]);
        outdeg[U[i]]++, indeg[V[i]]++;
    }
    vector<int> dp(N, 0);
    for (int v = 0; v < N; v++) {
        for (auto [to, w] : dag[v]) dp[to] = max(dp[to], dp[v] + w);
    }
    int D = dp[N-1];

    MinCostTension<long long, long long> opt(N);
    const long long INF = 1LL << 40;
    for (int v = 0; v < N; v++) opt.add_single_coef(v, outdeg[v] - indeg[v]);
    opt.add_tension_constraint(0, N-1, INF, D);
    for (int i = 0; i < M; i++) opt.add_tension_constraint(V[i], U[i], INF, -W[i]);
    auto [flag, cost] = opt.solve();
    auto res = -(sum + cost);  // 元々最大化問題だった
    auto ans = opt.reconstruct();
    cout << res << endl;
}

// AtCoder ABC 347 G - Grid Coloring 2
/*
    目的関数はテンションについての以下の区分線形凸関数
    ・最小値：0
    ・左側の変化点：{(0, 1), (-1, 2), (-2, 2), (-3, 2), (-4, 2)}
    ・右側の変化点：{(0, 1), (1, 2), (2, 2), (3, 2), (4, 2)}
    （傾きが 1, 3, 5, 7, 9 と変化していくので、その差分をとっている）

    制約条件
    ・p[v] - s <= 5
    ・s - p[v] <= -1
    ・値の決まっている p については p[v] - s <= (値), s - p[v] <= -(値)
*/
void ABC_347_G() {
    int N;
    cin >> N;
    const long long INF = 1LL << 30;
    vector<vector<int>> A(N, vector<int>(N));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) cin >> A[i][j];

    MinCostTension<long long, long long> opt(N * N + 1);
    int s = N * N;
    vector<pair<long long, long long>> left{{0, 1}, {-1, 2}, {-2, 2}, {-3, 2}, {-4, 2}};
    vector<pair<long long, long long>> right{{0, 1}, {1, 2}, {2, 2}, {3, 2}, {4, 2}};
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        if (i+1 < N) opt.add_tension_convex_function(i*N+j, (i+1)*N+j, 0, left, right);
        if (j+1 < N) opt.add_tension_convex_function(i*N+j, i*N+j+1, 0, left, right);
        if (A[i][j] == 0) {
            opt.add_tension_constraint(s, i*N+j, INF, 5);
            opt.add_tension_constraint(i*N+j, s, INF, -1);
        } else {
            opt.add_tension_constraint(s, i*N+j, INF, A[i][j]);
            opt.add_tension_constraint(i*N+j, s, INF, -A[i][j]);
        }
    }
    auto [flag, cost] = opt.solve();
    auto ans = opt.reconstruct();
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) cout << ans[i*N+j] - ans[s] << " ";
        cout << endl;
    }
}

//AtCoder ABC 397 G - Maximize Distance
/*
    最短距離を d 以上にできるか？　以下の最適化問題の答えが K 以下
    min
        Σ_{e} max(0, p[v] - p[u])
    s.t.
        p[N-1] - p[0] >= d
        p[v] - p[u] <= 1
*/
void ABC_397_G() {
    long long INF = 1LL << 30;
    long long N, M, K;
    cin >> N >> M >> K;
    vector<int> U(M), V(M), indeg(N, 0), outdeg(N, 0);
    for (int i = 0; i < M; i++) {
        cin >> U[i] >> V[i], U[i]--, V[i]--;
        indeg[V[i]]++, outdeg[U[i]]++;
    }

    long long low = -1, high = 1000;
    while (high - low > 1) {
        long long d = (high + low) / 2;
        MinCostTension<long long, long long> opt(N);
        opt.add_tension_constraint(N-1, 0, INF, -d);
        for (int i = 0; i < M; i++) {
            opt.add_tension_cost(U[i], V[i], 1, 0);
            opt.add_tension_constraint(U[i], V[i], INF, 1);
        }
        auto [flag, cost] = opt.solve();
        auto ans = opt.reconstruct();
        if (flag && cost <= K) low = d;
        else high = d;
    }
    cout << low << endl;
}

// HUPC 2020 H - Traditional Company (AOJ 3171)
/*
    T 以内にできるか？　以下を -X 以下にできるか？
    Σ_{i} -B[i](e[i]-s[i]) + (B[i]-A[i])C[i] + (A[i]-B[i])max(0, s[i]-e[i]+C[i])

    s[i] - s[i+1] <= 0 (for any i)
    e[i] - s[0] <= T (for any i)
    e[i] - s[i] >= 1 (for any i)

    naka[i] = 1 のとき
        e[U[i]] - s[V[i]] >= 1
    naka[i] = -1 のとき
        s[V[i]] >= e[U[i]]
*/
void AOJ_3171() {
    long long N, M, X, INF = 1LL<<40, LIM = 1LL << 30;
    cin >> N >> M >> X;
    vector<long long> A(N), B(N), C(N), U(M), V(M), naka(M);
    for (int i = 0; i < N; i++) cin >> A[i] >> B[i] >> C[i];
    for (int i = 0; i < M; i++) cin >> U[i] >> V[i] >> naka[i], U[i]--, V[i]--;

    long long low = -1, high = LIM;
    while (high - low > 1) {
        long long T = (low + high) / 2;
        MinCostTension<long long, long long> opt(N * 2);
        for (int i = 0; i < N; i++) {
            opt.add_cost((B[i] - A[i]) * C[i]);
            opt.add_single_coef(i, B[i]);
            opt.add_single_coef(i+N, -B[i]);
            opt.add_tension_cost(i+N, i, A[i]-B[i], -C[i]);
            if (i+1 < N) opt.add_tension_constraint(i+1, i, INF, 0);
            opt.add_tension_constraint(0, i+N, INF, T);
            opt.add_tension_constraint(i+N, i, INF, -1);
        }
        for (int i = 0; i < M; i++) {
            if (naka[i] == 1) opt.add_tension_constraint(U[i]+N, V[i], INF, -1);
            else opt.add_tension_constraint(V[i], U[i]+N, INF, 0);
        }
        auto [flag, cost] = opt.solve();
        if (flag && cost <= -X) high = T;
        else low = T;
    }
    cout << (high < LIM ? high : -1) << endl;
}

// AtCoder ABC 393 G - Unevenness
/*
    min: Σ_{u < v}(max(0, x[v] - x[u]) + max(0, x[u] - x[v]))
    s.t. Σ_{v}(max(0, x[v] - A[v]) + max(0, A[v] - x[v])) <= P/Q

    ラグランジュ緩和して
    max_{λ >= 0}:
        min: Σ_{u < v}(max(0, x[v] - x[u]) + max(0, x[u] - x[v]))
            + λ(Σ_{v}(max(0, x[v] - A[v]) + max(0, A[v] - x[v])) - P/Q)
*/
template<class T> struct SternBrocotTree {
    // binary search on Stern-Brocot Tree
    // return {l (= a/b), r (= c/d(} s.t. l: OK, r: NG
    // and a, b, c, d are maximized where a, b, c, d <= lim
    template<class Func> static tuple<T, T, T, T> binary_search(Func check, T lim) {
        assert(check(0, 1));
        assert(!check(1, 0));
        auto rec = [&](auto &&rec, bool which, T &a, T &b, T c, T d) -> void {
            if (a + c > lim || b + d > lim) return;
            if (check(a + c, b + d) == which) {
                a += c, b += d;
                rec(rec, which, a, b, c + c, d + d);
            }
            if (a + c <= lim && b + d <= lim && check(a + c, b + d) == which) a += c, b += d;
        };
        T a = 0, b = 1, c = 1, d = 0;
        while (a + c <= lim && b + d <= lim) {
            rec(rec, true, a, b, c, d);
            rec(rec, false, c, d, a, b);
        }
        return {a, b, c, d};
    }
};
template<class T = long long> struct frac {
    // gcd
    static T gcd(T a, T b) {
        a = max(a, -a), b = max(b, -b);
        while (b) {
            a %= b;
            swap(a, b);
        }
        return a;
    }

    // inner values
    T first, second;

    // constructor
    frac& normalize() {
        if (first == 0 && second != 0) {
            second = 1;
            return *this;
        }
        if (second == 0 && first != 0) {
            first = 1;
            return *this;
        }
        if (second < 0) first = -first, second = -second;
        T d = gcd(max(first, -first), second);
        if (d == 0) first = 0, second = 1;
        else first /= d, second /= d;
        return *this;
    }
    frac(const frac&) = default;
    frac& operator = (const frac&) = default;
    constexpr frac(T f = 0, T s = 1) : first(f), second(s) { 
        normalize(); 
    }
    constexpr frac& operator = (T a) { 
        *this = frac(a, 1); 
        return *this;
    }
    constexpr long double to_double() const {
        assert(second != 0);
        return (long double)(first) / (long double)(second);
    }
    friend constexpr long double to_double(const frac &r) {
        return r.to_double();
    }

    // comparison operators
    constexpr bool operator == (const frac &r) const {
        return this->first == r.first && this->second == r.second;
    }
    constexpr bool operator != (const frac &r) const {
        return this->first != r.first || this->second != r.second;
    }
    constexpr bool operator < (const frac &r) const {
        return this->first * r.second < this->second * r.first;
    }
    constexpr bool operator > (const frac &r) const {
        return this->first * r.second > this->second * r.first;
    }
    constexpr bool operator <= (const frac &r) const {
        return this->first * r.second <= this->second * r.first;
    }
    constexpr bool operator >= (const frac &r) const {
        return this->first * r.second >= this->second * r.first;
    }
    
    // arithmetic operators
    constexpr frac& operator += (const frac &r) {
        this->first = this->first * r.second + this->second * r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator -= (const frac &r) {
        this->first = this->first * r.second - this->second * r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator *= (const frac &r) {
        this->first *= r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator /= (const frac &r) {
        this->first *= r.second;
        this->second *= r.first;
        this->normalize();
        return *this;
    }
    constexpr frac operator + () const { return frac(*this); }
    constexpr frac operator - () const { return frac(0) - frac(*this); }
    constexpr frac operator + (const frac &r) const { return frac(*this) += r; }
    constexpr frac operator - (const frac &r) const { return frac(*this) -= r; }
    constexpr frac operator * (const frac &r) const { return frac(*this) *= r; }
    constexpr frac operator / (const frac &r) const { return frac(*this) /= r; }
    friend constexpr ostream& operator << (ostream &os, const frac<T> &x) {
        os << x.first; 
        if (x.second != 1) os << "/" << x.second;
        return os;
    }
};
void ABC_393_G() {
    using i128 = __int128_t;
    using FR = frac<i128>;
    using SBT = SternBrocotTree<long long>;
    long long N, P, Q;
    cin >> N >> P >> Q;
    FR K(P, Q);
    vector A(N, vector(N, 0LL));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) cin >> A[i][j];

    // 目的関数の値
    auto calc_obj = [&](const vector<FR> &x) -> FR {
        FR res = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            if (i+1 < N) res += max(x[i*N+j] - x[(i+1)*N+j], x[(i+1)*N+j] - x[i*N+j]);
            if (j+1 < N) res += max(x[i*N+j] - x[i*N+j+1], x[i*N+j+1] - x[i*N+j]);
        }
        return res;
    };

    // 小さい λ では負の値になり、ある程度大きい λ では 0 になる。0 になる瞬間が最適解。
    auto calc_penalty = [&](const vector<FR> &x) -> FR {
        FR res = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            res += max(x[i*N+j] - FR(A[i][j]), FR(A[i][j]) - x[i*N+j]);
        }
        return res - K;
    };

    // Stern-Brocot 木上の二分探索を実施する
    /*
    min: Σ_{u < v}(max(0, x[v] - x[u]) + max(0, x[u] - x[v]))
            + λ(Σ_{v}(max(0, x[v] - A[v]) + max(0, A[v] - x[v])) - P/Q)
    */
    auto optimize = [&](FR r, vector<FR> &x) -> FR {
        MinCostTension<FR, FR> opt(N * N + 1);
        int s = N * N;
        opt.add_cost(-r * K);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            if (i+1 < N) {
                int u = i*N+j, v = (i+1)*N+j;
                opt.add_tension_cost(u, v, FR(1), FR(0));
                opt.add_tension_cost(v, u, FR(1), FR(0));
            }
            if (j+1 < N) {
                int u = i*N+j, v = i*N+j+1;
                opt.add_tension_cost(u, v, FR(1), FR(0));
                opt.add_tension_cost(v, u, FR(1), FR(0));
            }
            int v = i*N+j;
            opt.add_tension_cost(s, v, r, FR(A[i][j]));
            opt.add_tension_cost(v, s, r, -FR(A[i][j]));
        }
        auto [flag, cost] = opt.solve();
        auto ans = opt.reconstruct();
        for (int i = 0; i < N * N; i++) x[i] = ans[i] - ans[s];
        return cost;
    };
    auto check = [&](i128 a, i128 b) -> bool {
        FR r(a, b);
        vector<FR> x(N*N);
        if (a == 1 && b == 0) return false;
        auto cost = optimize(r, x);
        return calc_penalty(x) > 0;
    };
    
    vector<FR> pre_x(N*N), nex_x(N*N), x(N*N);
    if (!check(0, 1)) {
        vector<FR> x(N*N);
        optimize(FR(0, 1), x);
    } else {
        auto [a, b, c, d] = SBT::binary_search(check, 10000000000000LL);
        FR pre_r(a, b), nex_r(c, d);
        auto pre_all_cost = optimize(pre_r, pre_x);
        auto nex_all_cost = optimize(nex_r, nex_x);
        auto pre_obj = calc_obj(pre_x), nex_obj = calc_obj(nex_x);
        auto pre_penalty = calc_penalty(pre_x), nex_penalty = calc_penalty(nex_x);
        for (int i = 0; i < N*N; i++) {
            x[i] = (pre_x[i]*(-nex_penalty) + nex_x[i]*pre_penalty) / (pre_penalty - nex_penalty);
        }
    }
    auto res = calc_obj(x);
    cout << fixed << setprecision(20) << to_double(res) << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) cout << to_double(x[i*N+j]) << " ";
        cout << endl;
    }
}


int main() {
    //AOJ_2230();
    //ABC_347_G();
    //ABC_397_G();
    //AOJ_3171();
    ABC_393_G();
}