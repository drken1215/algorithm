//
// 最小費用 b-flow by primal-dual (負閉路 NG)
//
// example:
//   ABC 421 G - Increase to make it Increasing (N <= 300)
//   (primal-dual: 4 ms, cost-scaling: 24 ms, simplex: 2 ms)
//     https://atcoder.jp/contests/abc421/tasks/abc421_g
//
//   KUPC 2014 I - Rain (N <= 10^4, ただし K <= 15 によって流量が小さい)
//   (primal-dual: 35 ms, cost-scaling: TLE, simplex: 681 ms)
//     https://atcoder.jp/contests/kupc2014/tasks/kupc2014_i
//
//   JAG 夏合宿 2013 Day4 I - Multi Path Story (N <= 1000)
//   (primal-dual: 18 ms, cost-scaling: 2590 ms, simplex: 6 ms)
//     https://onlinejudge.u-aizu.ac.jp/problems/2627
//
//   Educational Codeforces Round 80 F. Red-Blue Graph (N <= 200)
//   (primal-dual: 62 ms, cost-scaling: 46 ms, simplex: 62 ms)
//     https://codeforces.com/contest/1288/problem/F
//


#include <bits/stdc++.h>
using namespace std;


// output stream
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class S, class T> ostream& operator << (ostream &s, const pair<S, T> &P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 2> &P)
{ return s << '<' << P[0] << "," << P[1] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 3> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 4> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << "," << P[3] << '>'; }
template<class T> ostream& operator << (ostream &s, const vector<string> &P)
{ for (int i = 0; i < P.size(); ++i) { s << P[i] << endl; } return s; }
template<class T> ostream& operator << (ostream &s, const vector<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const deque<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const vector<vector<T>> &P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, const set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const multiset<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const unordered_set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const unordered_map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }


//--------------------------------//
// max flow, min-cost flow
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
    FLOW augment(int s, int t, vector<FlowEdge<FLOW>> &path, FLOW up_flow = numeric_limits<FLOW>::max()) {
        vector<bool> seen(size(), false);
        auto dfs = [&](auto &&dfs, int v, vector<FlowEdge<FLOW>> &path, FLOW up_flow) -> FLOW {
            if (v == t) return up_flow;
            seen[v] = true;
            for (int i = 0; i < (int)list[v].size(); i++) {
                FlowEdge<FLOW> &e = list[v][i], &re = get_rev_edge(e);
                if (seen[e.to] || e.cap <= 0) continue;
                FLOW flow = dfs(dfs, e.to, path, min(up_flow, e.cap));
                if (flow > 0) {
                    e.cap -= flow, e.flow += flow;
                    re.cap += flow, re.flow -= flow;
                    path.emplace_back(e);
                    return flow;
                }
            }  
            return FLOW(0); 
        };
        path.clear();
        FLOW res = dfs(dfs, s, path, up_flow);
        reverse(path.begin(), path.end());
        return res;
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
        assert(0 <= i && i < (int)list.size());
        return list[i];
    }
    const vector<FlowCostEdge<FLOW, COST>> &operator [] (int i) const {
        assert(0 <= i && i < (int)list.size());
        return list[i];
    }
    size_t size() const noexcept {
        return list.size();
    }
    size_t size_edegs() const noexcept {
        return pos.size();
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
        assert(0 <= from && from < (int)list.size() && 0 <= to && to < (int)list.size());
        assert(cap >= 0);
        int from_id = (int)list[from].size(), to_id = (int)list[to].size();
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, 0, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, 0, cap, -cost));
        if (cost < 0) include_negative_edge = true;
    }
    void add_edge(int from, int to, FLOW cap, FLOW rcap, COST cost) {
        assert(0 <= from && from < (int)list.size() && 0 <= to && to < (int)list.size());
        assert(cap >= 0);
        int from_id = (int)list[from].size(), to_id = (int)list[to].size();
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, rcap, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, rcap, cap, -cost));
        if (cost < 0) include_negative_edge = true;
    }
    void add_bidirected_edge(int from, int to, FLOW cap, COST cost) {
        assert(0 <= from && from < (int)list.size() && 0 <= to && to < (int)list.size());
        assert(cap >= 0);
        add_edge(from, to, cap, cap, cost);
    }

    // find initial potential (to resolve initial negative-edge)
    // pot[v] := potential (e.cost + pot[e.from] - pos[e.to] >= 0)
    bool calc_potential_dag() {
        pot.assign(size(), 0);
        vector<int> deg(size(), 0), st;
        for (int v = 0; v < (int)size(); v++) for (const auto &e : list[v]) deg[e.to] += (e.cap > 0);
        st.reserve(size());
        for (int v = 0; v < (int)size(); v++) if (!deg[v]) st.emplace_back(v);
        for (int i = 0; i < (int)size(); i++) {
            if ((int)st.size() == i) return false;  // not DAG
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
        for (int v = 0; v < (int)size(); v++) que.push(v), inque[v] = true;
        while (!que.empty()) {
            int cur = que.front();
            que.pop();
            inque[cur] = false;
            if (cnt[cur] > (int)size()) return false;  // include negative-cycle
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
            if (is_cycle) cycles.push_back(std::move(seq));
            else paths.push_back(std::move(seq));
        };

        // Phase 1: extract all cycles and make graph DAG
        const int NOTSEEN = 0, INSTACK = 1, FINISH = 2;
        vector<int> color(list.size(), NOTSEEN);
        vector<int> pos_in_stack(list.size(), -1);
        vector<pair<int, int>> stk;
        auto dfs = [&](auto &&dfs, int v) -> bool {
            color[v] = INSTACK;
            pos_in_stack[v] = (int)stk.size();
            for (int i = 0; i < (int)fg[v].size(); i++) {
                if (fg[v][i].rem <= 0) continue;
                int u = fg[v][i].to;
                if (color[u] == INSTACK) {
                    vector<pair<int,int>> route;
                    for (int k = pos_in_stack[u]; k < (int)stk.size(); k++) {
                        route.push_back(stk[k]);
                    }
                    route.push_back({v, i});
                    build(route, true);
                    return true;
                }
                if (color[u] == NOTSEEN) {
                    stk.push_back({v, i});
                    if (dfs(dfs, u)) return true;
                    stk.pop_back();
                }
            }
            color[v] = FINISH;
            pos_in_stack[v] = -1;
            return false;
        };
        while (true) {
            fill(color.begin(), color.end(), NOTSEEN);
            stk.clear();
            bool found = false;
            for (int v = 0; v < (int)list.size() && !found; v++) {
                if (color[v] == NOTSEEN && dfs(dfs, v)) found = true;
            }
            if (!found) break;
        }

        // Phase 2: find all s-t paths
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

// min-cost max-flow (<= limit_flow), slope ver.
template<class FLOW, class COST> vector<pair<FLOW, COST>>
MinCostFlowSlope(FlowCostGraph<FLOW, COST> &G, int S, int T, FLOW limit_flow)
{
    // result values
    FLOW cur_flow = 0;
    COST cur_cost = 0, pre_cost = numeric_limits<COST>::max() / 2;
    vector<pair<FLOW, COST>> res;
    res.emplace_back(cur_flow, cur_cost);
    
    // intermediate values
    vector<COST> dist((int)G.size(), numeric_limits<COST>::max() / 2);
    vector<int> prevv((int)G.size(), -1), preve((int)G.size(), -1);
    
    // dual
    auto dual_step = [&]() -> bool {
        dist.assign((int)G.size(), numeric_limits<COST>::max() / 2);
        dist[S] = 0;
        priority_queue<pair<COST,int>, vector<pair<COST,int>>, greater<pair<COST,int>>> que;
        que.emplace(0, S);
        while (!que.empty()) {
            auto [cur, v] = que.top();
            que.pop();
            if (cur > dist[v]) continue;
            for (int i = 0; i < (int)G[v].size(); i++) {
                const auto &e = G[v][i];
                COST add = e.cost + G.pot[v] - G.pot[e.to];
                if (e.cap > 0 && dist[e.to] > dist[v] + add) {
                    dist[e.to] = dist[v] + add;
                    prevv[e.to] = v;
                    preve[e.to] = i;
                    que.emplace(dist[e.to], e.to);
                }
            }
        }
        return dist[T] < numeric_limits<COST>::max() / 2;
    };
    
    // primal
    auto primal_step = [&]() -> void {
        for (int v = 0; v < G.size(); v++) {
            if (dist[v] < numeric_limits<COST>::max() / 2) G.pot[v] += dist[v];
            else G.pot[v] = numeric_limits<COST>::max() / 2;
        }
        FLOW flow = limit_flow - cur_flow;
        COST cost = G.pot[T] - G.pot[S];
        for (int v = T; v != S; v = prevv[v]) {
            flow = min(flow, G[prevv[v]][preve[v]].cap);
        }
        for (int v = T; v != S; v = prevv[v]) {
            FlowCostEdge<FLOW, COST> &e = G[prevv[v]][preve[v]];
            FlowCostEdge<FLOW, COST> &re = G.get_rev_edge(e);
            e.cap -= flow, e.flow += flow;
            re.cap += flow, re.flow -= flow;
        }
        cur_flow += flow;
        cur_cost += flow * cost;
        if (pre_cost == cost) res.pop_back();
        res.emplace_back(cur_flow, cur_cost);
        pre_cost = cost;
    };

    // initialize potential
    assert(G.init_potential());
    
    // primal-dual
    while (cur_flow < limit_flow) {
        if (!dual_step()) break;
        primal_step();
    }
    return res;
}

// min-cost max-flow, slope ver.
template<class FLOW, class COST> vector<pair<FLOW, COST>>
MinCostFlowSlope(FlowCostGraph<FLOW, COST> &G, int S, int T)
{
    return MinCostFlowSlope(G, S, T, numeric_limits<FLOW>::max());
}

// min-cost max-flow (<= limit_flow)
template<class FLOW, class COST> pair<FLOW, COST>
MinCostFlow(FlowCostGraph<FLOW, COST> &G, int S, int T, FLOW limit_flow)
{
    return MinCostFlowSlope(G, S, T, limit_flow).back();
}

// min-cost max-flow (<= limit_flow)
template<class FLOW, class COST> pair<FLOW, COST>
MinCostFlow(FlowCostGraph<FLOW, COST> &G, int S, int T)
{
    return MinCostFlow(G, S, T, numeric_limits<FLOW>::max());
}


//--------------------------------//
// b-flow
//--------------------------------//

// Minimum Cost b-flow (by primal-dual, negative cycle is NG)
template<class FLOW, class COST> struct MinCostBFlowByPrimalDual {
    // inner values
    int N;
    FlowCostGraph<FLOW, COST> G;
    vector<FLOW> dss;  // demand (< 0) and supply (> 0)
    vector<COST> dual;

    // constructor
    explicit MinCostBFlowByPrimalDual(int n) : N(n), G(n + 2), dss(n, 0) {}

    // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        G.add_edge(from, to, cap, cost);
    }
    void set_ds(int v, FLOW ds) {
        assert(0 <= v && v < N);
        dss[v] = ds;
    }
    void set_ds(const vector<FLOW> &vds) {
        assert((int)vds.size() == N);
        dss = vds;
    }

    // getter
    FlowCostEdge<FLOW, COST> &get_edge(int i) {
        return G.get_edge(i);
    }
    const FlowCostEdge<FLOW, COST> &get_edge(int i) const {
        return G.get_edge(i);
    }
    vector<FlowCostEdge<FLOW, COST>> get_edges() const {
        return G.get_edges();
    }
    COST get_dual(int v) const {
        return dual[v];
    }
    vector<COST> get_duals() const {
        return dual;
    }

    // solver
    pair<bool, COST> solve(bool calc_potential = false) {
        // dss treatment
        int s = N, t = N + 1;
        FLOW ssum = 0, tsum = 0;
        for (int v = 0; v < N; v++) {
            if (dss[v] > 0) ssum += dss[v], G.add_edge(s, v, dss[v], COST(0));
            else if (dss[v] < 0) tsum -= dss[v], G.add_edge(v, t, -dss[v], COST(0));
        }

        // feasibility check
        if (ssum != tsum) return {false, COST(0)};
        
        // min-cost flow
        auto [maxflow, mincost] = MinCostFlow(G, s, t, ssum);
        if (maxflow < ssum) return {false, COST(0)};

        // find dual
        if (calc_potential) {
            G.calc_potential();
            dual = G.pot;
            dual.pop_back(), dual.pop_back();  // eliminate s, t
        }
        return {true, mincost};
    }
};

// b-flow manager
template<class FLOW, class COST> struct MinCostBFlow {
    // Edge
    struct InnerEdge {
        int from, to;
        FLOW lower_cap, upper_cap, flow;
        COST cost;
        InnerEdge(int from_, int to_, FLOW lower_, FLOW upper_, COST cost_)
            : from(from_), to(to_), lower_cap(lower_), upper_cap(upper_), flow(0), cost(cost_) {}
        friend ostream& operator << (ostream& s, const InnerEdge& e) {
            return s << e.from << "->" << e.to 
            << " (" << e.flow << "/" << e.lower_cap << "~" << e.upper_cap << ", " << e.cost << ")";
        }
    };

    // inner values
    int N;
    vector<InnerEdge> edges;
    vector<FLOW> lower_dss, upper_dss, dss;  // demand (< 0) and supply (> 0)
    vector<COST> dual;
    
    // constructor
    explicit MinCostBFlow(int n = 0) : N(n), lower_dss(n, 0), upper_dss(n, 0), dss(n, 0) {}

    // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        edges.push_back(InnerEdge(from, to, 0, cap, cost));
    }
    void add_edge(int from, int to, FLOW lower_cap, FLOW upper_cap, COST cost) {
        assert(lower_cap <= upper_cap);
        edges.push_back(InnerEdge(from, to, lower_cap, upper_cap, cost));
    }
    void set_ds(int v, FLOW ds) {
        assert(0 <= v && v < N);
        lower_dss[v] = ds, upper_dss[v] = ds;
    }
    void set_ds(int v, FLOW lower_ds, FLOW upper_ds) {
        assert(0 <= v && v < N);
        assert(lower_ds <= upper_ds);
        lower_dss[v] = lower_ds, upper_dss[v] = upper_ds;
    }

    // getter
    InnerEdge &get_edge(int i) {
        return edges[i];
    }
    const InnerEdge &get_edge(int i) const {
        return edges[i];
    }
    vector<InnerEdge> get_edges() const {
        return edges;
    }
    COST get_dual(int v) const {
        return dual[v];
    }
    vector<COST> get_duals() const {
        return dual;
    }

    // solver
    bool pre_compute() {
        bool need_super_node = false;
        for (int v = 0; v < N; v++) {
            if (lower_dss[v] == upper_dss[v]) dss[v] = lower_dss[v];
            else need_super_node = true;
        }

        // lower_ds, upper_ds -> strict ds
        if (need_super_node) {
            int super = N;
            dss.assign(N + 1, 0);
            for (int v = 0; v < N; v++) {
                if (lower_dss[v] >= 0) {
                    add_edge(super, v, lower_dss[v], upper_dss[v], 0);
                } else if (upper_dss[v] < 0) {
                    add_edge(v, super, -upper_dss[v], -lower_dss[v], 0);
                } else {
                    add_edge(super, v, upper_dss[v], 0);
                    add_edge(v, super, -lower_dss[v], 0);
                }
            }
        }

        // push lower_cap
        for (const auto &e : edges) {
            dss[e.to] += e.lower_cap, dss[e.from] -= e.lower_cap;
        }
        return need_super_node;
    }
    pair<bool, COST> solve(const string solver = "network_simplex", bool calc_potential = false) {
        bool need_super_node = pre_compute();
        COST res = 0;
        if (solver == "primal_dual") {
            MinCostBFlowByPrimalDual<FLOW, COST> G(N + (int)need_super_node);
            G.set_ds(dss);
            for (const auto &e : edges) G.add_edge(e.from, e.to, e.upper_cap - e.lower_cap, e.cost);
            auto [feasible, mincost] = G.solve(calc_potential);
            if (!feasible) return {false, COST(0)};
            for (int i = 0; i < (int)edges.size(); i++) {
                auto &e = edges[i];
                const auto &ge = G.get_edge(i);
                e.flow = e.upper_cap - ge.cap;
                res += e.flow * e.cost;
            }
            if (calc_potential) {
                dual = G.get_duals();
                if (need_super_node) dual.pop_back();
            }
        }
        return {true, res};
    }
};


//------------------------------//
// Solver
//------------------------------//

// ABC 421 G - Increase to make it Increasing
void ABC_421_G(const string solver) {
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
    auto [flag, cost] = G.solve(solver);
    cout << (flag ? cost : -1) << endl;
}

// KUPC 2014 I - Rain
void KUPC_2014_I(const string solver) {
    int N, M, K, INF = 1LL<<20;
    cin >> N >> M >> K;
    vector<long long> C(K), A(M), B(M), D(M), b(N, 0); 
    for (int i = 0; i < K; i++) cin >> C[i], C[i]--;
    for (int i = 0; i < M; i++) cin >> A[i] >> B[i] >> D[i], A[i]--, B[i]--;
    for (int i = 0; i < K; i++) b[A[C[i]]]++, b[B[C[i]]]--;
    MinCostBFlow<int, long long> G(N);
    for (int i = 0; i < N; i++) G.set_ds(i, b[i]);
    for (int i = 0; i < M; i++) G.add_edge(B[i], A[i], INF, D[i]);
    auto [flag, cost] = G.solve(solver);
    cout << (flag ? cost : -1) << endl;
}

// JAG 夏合宿 2013 Day4 I - Multi Path Story (AOJ 2627)
void AOJ_2627(const string solver) {
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
    auto [flag, cost] = G.solve(solver);
    cout << cost << endl;
}

// Educational Codeforces Round 80 F. Red-Blue Graph
void EducationalCodeforces80_F(const string solver) {
    long long L, R, M, costR, costB, INF = 10000;
    string sl, sr;
    cin >> L >> R >> M >> costR >> costB >> sl >> sr;
    MinCostBFlow<long long, long long> G(L + R);
    for (int i = 0; i < L; i++) {
        if (sl[i] == 'R') G.set_ds(i, 1, INF);
        else if (sl[i] == 'B') G.set_ds(i, -INF, -1);
        else G.set_ds(i, -INF, INF);
    }
    for (int j = 0; j < R; j++) {
        if (sr[j] == 'R') G.set_ds(j+L, -INF, -1);
        else if (sr[j] == 'B') G.set_ds(j+L, 1, INF);
        else G.set_ds(j+L, -INF, INF);
    }
    for (int i = 0; i < M; i++) {
        int u, v;
        cin >> u >> v, u--, v--;
        G.add_edge(u, v+L, 1, costR);
        G.add_edge(v+L, u, 1, costB);
    }
    auto [flag, mincost] = G.solve(solver);
     if (!flag) {
        cout << -1 << endl;
        return;
     }
    auto es = G.get_edges();
    string res(M, 'U');
    for (int i = 0; i < M; i++) {
        auto e = es[i*2], re = es[i*2+1];
        if (e.flow == 1) res[i] = 'R';
        else if (re.flow == 1) res[i] = 'B';
    }
    cout << mincost << endl;
    cout << res << endl;
}


int main() {
    ABC_421_G("primal_dual");
    //KUPC_2014_I("primal_dual");
    //AOJ_2627("primal_dual");
    //EducationalCodeforces80_F("primal_dual");
}