//
// 最小費用 b-flow by three methods:
//    ・primal-dual (負閉路 NG)
//    ・cost-scaling
//    ・network simplex (多くの場合最速)
//
// verified
//   Yosupo Library Checker - Minimum Cost b-flow (N <= 100)
//   (primal-dual: 負閉路 NG, cost-scaling: 15 ms, simplex: 2 ms)
//     https://judge.yosupo.jp/problem/min_cost_b_flow
//
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
//   UTPC 2011 H - キャッシュ戦略 (N <= 10^4, 本来は流量 M <= 10 だがそれを活かさない解法をしている）
//   (primal-dual: TLE, cost-scaling: TLE, simplex: 1130 ms)
//   #pragma GCC optimize("Ofast") を入れるとむしろ遅くなる！！！
//     https://atcoder.jp/contests/utpc2011/tasks/utpc2011_8
//
//   Codeforces Round 826 (Div. 3) G. Kirill and Company 
//   (primal-dual: 負閉路 NG, cost-scaling: 1046 ms, simplex: 62 ms)
//     https://codeforces.com/contest/1741/problem/G
//
//   AtCoder ABC 393 G - Unevenness (for using frac<i128> and dual)
//   (primal-dual: 負閉路 NG, cost-scaling: 830 ms, simplex: 116 ms)
//     https://atcoder.jp/contests/abc393/tasks/abc393_g
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
    vector<COST> pot; // pot[v] := potential (e.cost + pot[e.from] - pos[e.to] >= 0)
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

// Min Cost Circulation Flow by Cost-Scaling 
template<class FLOW, class COST> COST MinCostCirculation(FlowCostGraph<FLOW, COST> &G) {
    const int N = (int)G.size();
    const COST SCALE = N + 1;
    COST eps = 1;
    vector<FLOW> balance(G.size(), 0);
    vector<COST> price(G.size(), 0);
    
    auto reduced_cost = [&](const FlowCostEdge<FLOW, COST> &e) -> COST {
        return e.cost * SCALE - price[e.from] + price[e.to];
    };

    auto ConstructGaux = [&]() -> void {
        vector<bool> visited(G.size(), false);
        vector<int> st;
        st.reserve(N);
        for (int s = 0; s < N; s++) {
            if (balance[s] <= 0 || visited[s]) continue;
            visited[s] = true;
            st.push_back(s);
            while (!st.empty()) {
                int v = st.back();
                st.pop_back();
                for (const auto &e : G[v]) {
                    if (e.cap <= 0 || reduced_cost(e) >= 0 || visited[e.to]) continue;
                    visited[e.to] = true;
                    st.push_back(e.to);
                }
            }
        }
        for (int v = 0; v < G.size(); ++v) if (visited[v]) price[v] += eps;
    };

    auto augment_blocking_flow = [&]() -> bool {
        vector<int> iter(N, 0);
        auto augment = [&](auto &&augment, int v, FLOW flow) -> FLOW {
            if (balance[v] < 0) {
                FLOW dif = min(flow, -balance[v]);
                balance[v] += dif;
                return dif;
            }
            for (int &i = iter[v]; i < (int)G[v].size(); i++) {
                auto &e = G[v][i];
                if (e.cap <= 0 || reduced_cost(e) >= 0) continue;
                FLOW dif = augment(augment, e.to, min(flow, e.cap));
                if (dif <= 0) continue;
                auto &re = G.get_rev_edge(e);
                e.cap -= dif, e.flow += dif;
                re.cap += dif, re.flow -= dif;
                return dif;
            }
            return FLOW(0);
        };
        bool finish = true;
        for (int v = 0; v < N; ++v) {
            while (balance[v] > 0) {
                FLOW f = augment(augment, v, balance[v]);
                if (f <= 0) break;
                balance[v] -= f;
            }
            if (balance[v] > 0) finish = false;
        }
        return finish;
    };

    // eps init
    COST need = 0;
    for (int v = 0; v < N; v++) {
        for (const auto &e : G[v]) {
            if (e.cap <= 0) continue;
            need = max(need, -e.cost * SCALE);
        }
    }
    while (eps < need) eps *= 2;

    // cost scaling
    while (eps > 1) {
        eps /= 2;
        for (int v = 0; v < N; v++) {
            for (int i = 0; i < (int)G[v].size(); i++) {
                auto &e = G[v][i];
                if (e.cap <= 0 || reduced_cost(e) >= 0) continue;
                auto &re = G.get_rev_edge(e);
                FLOW f = e.cap;
                balance[e.from] -= f, balance[e.to] += f;
                e.cap -= f, e.flow += f;
                re.cap += f, re.flow -= f;
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

// Minimum Cost b-flow (by cost-scaling min-cost circulation)
template<class FLOW, class COST> struct MinCostBFlowByCostScaling {
    // inner Edge
    struct InnerEdge {
        int from, to;
        FLOW cap;
        COST cost;
        InnerEdge(int from_, int to_, FLOW cap_, COST cost_) : from(from_), to(to_), cap(cap_), cost(cost_) {}
        friend ostream& operator << (ostream& s, const InnerEdge& e) {
            return s << e.from << " -> " << e.to << " (" << e.cap << ", " << e.cost << ")";
        }
    };

    // inner values
    int N;
    FlowCostGraph<FLOW, COST> G;
    vector<InnerEdge> edges;
    vector<FLOW> dss;  // demand (< 0) and supply (> 0)
    vector<COST> dual;

    // constructor
    explicit MinCostBFlowByCostScaling(int n = 0) : N(n), G(n), dss(n, 0) {}

    // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        edges.push_back(InnerEdge(from, to, cap, cost));
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
    pair<bool, COST> solve(bool calc_potential = true) {
        // push s-t flow
        FlowGraph<FLOW> preG(N + 2);
        int s = N, t = N + 1;
        for (const auto &e : edges) preG.add_edge(e.from, e.to, e.cap);
        FLOW ssum = 0, tsum = 0;
        for (int v = 0; v < N; v++) {
            if (dss[v] > 0) ssum += dss[v], preG.add_edge(s, v, dss[v]);
            else if (dss[v] < 0) tsum -= dss[v], preG.add_edge(v, t, -dss[v]);
        }

        // feasibility check
        if (ssum != tsum) return {false, COST(0)};
        if (Dinic(preG, s, t) < ssum) return {false, COST(0)};

        // come down to min-cost circulation
        for (int i = 0; i < (int)edges.size(); i++) {
            const auto &e = edges[i];
            const auto &ge = preG.get_edge(i);
            G.add_edge(ge.from, ge.to, ge.cap, ge.flow, e.cost);
        }
        COST mincost = MinCostCirculation(G);

        // find dual
        if (calc_potential) {
            G.calc_potential();
            dual = G.pot;
        }
        return {true, mincost};
    }
};

// Network Simplex Method
template<class FLOW, class COST> struct MinCostBFlowByNetworkSimplex {
    // inner Edge
    struct InnerEdge {
        int from, to;
        FLOW cap;
        COST cost;
        InnerEdge(int from_, int to_, FLOW cap_, COST cost_) : from(from_), to(to_), cap(cap_), cost(cost_) {}
    };
    struct Parent {
        int p, e;
        FLOW up, down;
    };

    // inner values
    int N, original_edge_size;
    vector<InnerEdge> edges;
    vector<FLOW> dss;  // demand (< 0) and supply (> 0)
    bool feasible;
    COST total_cost;
    vector<COST> dual;

    // intermediate results
    int BUCKET_SIZE, MINOR_LIMIT;
    vector<Parent> parents;
    vector<int> depth, nex, pre, candidates;

    // constructor
    explicit MinCostBFlowByNetworkSimplex(int n = 0) : N(n), dss(n) {}

    // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        edges.emplace_back(from, to, cap, cost);
        edges.emplace_back(to, from, 0, -cost);
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
    FLOW get_flow(int i) const {
        return edges[(i * 2) ^ 1].cap;
    }
    COST get_dual(int v) const {
        return dual[v];
    }
    vector<COST> get_duals() const {
        return dual;
    }

    // solver
    pair<bool, COST> solve() {
        BUCKET_SIZE = max(int(sqrt(double(edges.size())) * 0.2), 10);
        MINOR_LIMIT = max(int(BUCKET_SIZE * 0.1), 3);
        precompute();
        candidates.reserve(BUCKET_SIZE);
        int ei = 0;
        while (true) {
            for (int i = 0; i < MINOR_LIMIT; i++) if (!minor()) break;
            COST best = 0;
            int best_ei = -1;
            candidates.clear();
            for (int i = 0; i < (int)edges.size(); i++) {
                if (edges[ei].cap > 0) {
                    COST clen = edges[ei].cost + dual[edges[ei ^ 1].to] - dual[edges[ei].to];
                    if (clen < 0) {
                        if (clen < best) best = clen, best_ei = ei;
                        candidates.push_back(ei);
                        if ((int)candidates.size() == BUCKET_SIZE) break;
                    }
                }
                ei++;
                if (ei == (int)edges.size()) ei = 0;
            }
            if (candidates.empty()) break;
            push_flow(best_ei);
        }
        if (!postcompute()) return {false, COST(-1)};
        else return {true, total_cost};
    }

    void connect(int a, int b) {
        nex[a] = b, pre[b] = a;
    }

    void precompute() {
        original_edge_size = (int)edges.size();
        dual.assign(N + 1, 0); 
        parents.resize(N), depth.assign(N + 1, 1); 
        nex.assign((N + 1) * 2, 0), pre.assign((N + 1) * 2, 0);
        COST inf_cost = 1;
        for (int i = 0; i < (int)edges.size(); i += 2) {
            inf_cost += (edges[i].cost >= 0 ? edges[i].cost : -edges[i].cost);
        }
        edges.reserve((int)edges.size() + N * 2);
        for (int i = 0; i < N; i++) {
            if (dss[i] >= 0) {
                edges.push_back(InnerEdge(i, N, 0, inf_cost));
                edges.push_back(InnerEdge(N, i, dss[i], -inf_cost));
                dual[i] = -inf_cost;
            } else {
                edges.push_back(InnerEdge(i, N, -dss[i], -inf_cost));
                edges.push_back(InnerEdge(N, i, 0, inf_cost));
                dual[i] = inf_cost;
            }
            int e = (int)edges.size() - 2;
            parents[i] = {N, e, edges[e].cap, edges[e ^ 1].cap};
        }
        depth[N] = 0;
        for (int i = 0; i < N + 1; i++) connect(i * 2, i * 2 + 1);
        for (int i = 0; i < N; i++) connect(i * 2 + 1, nex[N * 2]), connect(N * 2, i * 2);
    }

    bool postcompute() {
        for (int i = 0; i < N; i++) {
            edges[parents[i].e].cap = parents[i].up;
            edges[parents[i].e ^ 1].cap = parents[i].down;
        }
        feasible = true;
        for (int i = 0; i < N; i++) {
            int e = original_edge_size + i * 2;
            if (dss[i] >= 0) {
                if (edges[e ^ 1].cap > 0) feasible = false;
            } else {
                if (edges[e].cap > 0) feasible = false;
            }
        }
        if (!feasible) return false;
        total_cost = 0;
        for (int i = 0; i < (int)edges.size(); i += 2) {
            total_cost += edges[i ^ 1].cap * edges[i].cost;
        }
        dual.pop_back();
        return true;
    }

    void push_flow(int ei0) {
        int u0 = edges[ei0 ^ 1].to, v0 = edges[ei0].to, del_u = v0;
        FLOW f = edges[ei0].cap;
        COST clen = edges[ei0].cost + dual[u0] - dual[v0];
        bool del_u_side = true;
        int lca = get_lca(u0, v0, f, del_u_side, del_u);
        if (f > 0) {
            int u = u0, v = v0;
            while (u != lca) parents[u].up += f, parents[u].down -= f, u = parents[u].p;
            while (v != lca) parents[v].up -= f, parents[v].down += f, v = parents[v].p;
        }
        int u = u0, par = v0;
        auto p_caps = make_pair(edges[ei0].cap - f, edges[ei0 ^ 1].cap + f);
        COST p_diff = -clen;
        if (!del_u_side) {
            swap(u, par); 
            swap(p_caps.first, p_caps.second);
            p_diff *= -1;
        }
        int par_e = ei0 ^ (del_u_side ? 0 : 1);
        while (par != del_u) {
            int d = depth[par], idx = u * 2;
            while (idx != u * 2 + 1) {
                if (idx % 2 == 0) d++, dual[idx / 2] += p_diff, depth[idx / 2] = d;
                else d--;
                idx = nex[idx];
            }
            connect(pre[u * 2], nex[u * 2 + 1]);
            connect(u * 2 + 1, nex[par * 2]);
            connect(par * 2, u * 2);
            swap(parents[u].e, par_e);
            par_e ^= 1;
            swap(parents[u].up, p_caps.first); 
            swap(parents[u].down, p_caps.second);
            swap(p_caps.first, p_caps.second);
            int next_u = parents[u].p; 
            parents[u].p = par;
            par = u;
            u = next_u;
        }
        edges[par_e].cap = p_caps.first;
        edges[par_e ^ 1].cap = p_caps.second;
    }

    bool minor() {
        if (candidates.empty()) return false;
        COST best = 0;
        int best_ei = -1;
        int i = 0;
        while (i < int(candidates.size())) {
            int ei = candidates[i];
            if (edges[ei].cap <= 0) {
                swap(candidates[i], candidates.back());
                candidates.pop_back();
                continue;
            }
            COST clen = edges[ei].cost + dual[edges[ei ^ 1].to] - dual[edges[ei].to];
            if (clen >= 0) {
                swap(candidates[i], candidates.back());
                candidates.pop_back();
                continue;
            }
            if (clen < best) best = clen, best_ei = ei;
            i++;
        }
        if (best_ei == -1) return false;
        push_flow(best_ei);
        return true;
    }

    int get_lca(int u, int v, FLOW &flow, bool &del_u_side, int &del_u) {
        auto up_u = [&]() {
            if (parents[u].down < flow) flow = parents[u].down, del_u = u, del_u_side = true;
            u = parents[u].p;
        };
        auto up_v = [&]() {
            if (parents[v].up <= flow) flow = parents[v].up, del_u = v, del_u_side = false;
            v = parents[v].p;
        };
        if (depth[u] >= depth[v]) {
            int num = depth[u] - depth[v];
            for (int i = 0; i < num; i++) up_u();
        } else {
            int num = depth[v] - depth[u];
            for (int i = 0; i < num; i++) up_v();
        }
        while (u != v) up_u(), up_v();
        return u;
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
        } else if (solver == "cost_scaling") {
            MinCostBFlowByCostScaling<FLOW, COST> G(N + (int)need_super_node);
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
        } else if (solver == "network_simplex") {
            MinCostBFlowByNetworkSimplex<FLOW, COST> G(N + (int)need_super_node);
            G.set_ds(dss);
            for (const auto &e : edges) G.add_edge(e.from, e.to, e.upper_cap - e.lower_cap, e.cost);
            auto [feasible, mincost] = G.solve();
            if (!feasible) return {false, COST(0)};
            for (int i = 0; i < (int)edges.size(); i++) {
                auto &e = edges[i];
                e.flow = e.lower_cap + G.get_flow(i);
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

// Yosupo Libray Checker - Minimum Cost b-flow
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
void Yosupo_Minimum_Cost_b_flow(const string solver) {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);

    int N, M;
    cin >> N >> M;
    MinCostBFlow<long long, i128> G(N);
    vector<i128> B(N);
    for (int i = 0; i < N; i++) cin >> B[i], G.set_ds(i, B[i]);
    vector<int> s(M), t(M);
    vector<i128> l(M), u(M), c(M);
    for (int i = 0; i < M; i++) {
        cin >> s[i] >> t[i] >> l[i] >> u[i] >> c[i];
        G.add_edge(s[i], t[i], l[i], u[i], c[i]);
    }
    auto [exist, res] = G.solve(solver, true);
    if (!exist) cout << "infeasible" << '\n';
    else {
        const auto &dual = G.get_duals();
        const auto &es = G.get_edges();
        cout << res << '\n';
        for (auto v : dual) cout << v << '\n';
        for (const auto &e : es) cout << e.flow << '\n';
    }
}

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

// UTPC 2011 H - キャッシュ戦略
void UTPC2011_H(const string solver) {
    int M, N, K;
    cin >> M >> N >> K;
    vector<int> W(N), A(K);
    for (int i = 0; i < N; i++) cin >> W[i];
    for (int i = 0; i < K; i++) cin >> A[i], A[i]--;
    int s = K * 3, t = s + 1;
    MinCostBFlow<int, int> G(K * 3 + 2);
    for (int v = 0; v < K; v++) {
        G.add_edge(s, v, 1, W[A[v]]);
        G.add_edge(v+K*2, v, 1, W[A[v]]);  // 空の状態からは常に行ける
        G.add_edge(v+K, t, 1, 0);
        G.add_edge(v, v+K, 1, 1, 0);  // 流量下限も 1

        // 次の同じ色のボールが来るまでキープする場合
        for (int v2 = v+1; v2 < K; v2++) {
            if (A[v] == A[v2]) {
                G.add_edge(v+K, v2, 1, 0);
                break;
            }
        }

        // 箱からボールを取り出す頂点との絡み
        if (v+1 < K) {
            G.add_edge(v+K, (v+1)+K*2, 1, 0);
            G.add_edge(v+K*2, (v+1)+K*2, M, 0);
        }
    }
    G.add_edge(t, s, M, 0);
    auto [flag, mincost] = G.solve(solver);
    cout << mincost << endl;
}

// Codeforces Round 826 (Div. 3) G. Kirill and Company
void Codeforces826_G(const string solver) {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int T;
    cin >> T;
    while (T--) {
        int N, M, F, K, INF = 6;
        cin >> N >> M;
        vector<vector<int>> G(N);
        for (int i = 0; i < M; i++) {
            int a, b; cin >> a >> b; a--, b--;
            G[a].emplace_back(b), G[b].emplace_back(a);
        }
        vector<vector<int>> prev(N);
        vector<int> dp(N, -1);
        queue<int> que;
        dp[0] = 0;
        que.push(0);
        while (!que.empty()) {
            auto v = que.front(); que.pop();
            for (auto v2 : G[v]) {
                if (dp[v2] == -1) {
                    dp[v2] = dp[v] + 1;
                    prev[v2].emplace_back(v);
                    que.push(v2);
                } else if (dp[v2] == dp[v] + 1) {
                    prev[v2].emplace_back(v);
                }
            }
        }
        cin >> F;
        vector<int> where(F), allnum(N, 0), nothavenum(N, 0), havenum(N, 0);
        for (int i = 0; i < F; i++) {
            cin >> where[i], where[i]--;
            allnum[where[i]]++;
        }
        cin >> K;
        for (int i = 0; i < K; i++) {
            int v; cin >> v, v--;
            nothavenum[where[v]]++;
        }
        for (int v = 0; v < N; v++) havenum[v] = allnum[v] - nothavenum[v];

        int s = N * 2, t = N;
        MinCostBFlow<int, int> FG(N * 2 + 1);
        for (int v = 0; v < N; v++) {
            if (havenum[v] > 0) FG.add_edge(s, v, havenum[v], 0);
            if (nothavenum[v] > 0) FG.add_edge(v, v+N, 1, -nothavenum[v]);
            FG.add_edge(v, v+N, INF, 0);
            for (auto v2 : prev[v]) FG.add_edge(v+N, v2, INF, 0);
        }
        FG.add_edge(t, s, INF, 0);
        auto [flag, mincost] = FG.solve(solver);
        int res = K + mincost;
        cout << res << endl;
    }
}

// AtCoder ABC 393 G - Unevenness
template<class T> struct SternBrocotTree {
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
void ABC_393_G(const string solver) {
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
        MinCostBFlow<FR, FR> G(N * N + 1);
        int s = N * N;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            if (i+1 < N) {
                int u = i*N+j, v = (i+1)*N+j;
                G.add_edge(u, v, FR(1), FR(0));
                G.add_edge(v, u, FR(1), FR(0));
            }
            if (j+1 < N) {
                int u = i*N+j, v = i*N+j+1;
                G.add_edge(u, v, FR(1), FR(0));
                G.add_edge(v, u, FR(1), FR(0));
            }
            int v = i*N+j;
            G.add_edge(s, v, r, FR(A[i][j]));
            G.add_edge(v, s, r, -FR(A[i][j]));
        }
        auto [flag, cost] = G.solve(solver, true);  // dual も求める
        auto ans = G.dual;
        assert(ans.size() == N * N + 1);
        for (int i = 0; i < N * N; i++) x[i] = ans[i] - ans[s];
        return cost - r * K;
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
    //Yosupo_Minimum_Cost_b_flow("primal_dual");
    //Yosupo_Minimum_Cost_b_flow("cost_scaling");
    //Yosupo_Minimum_Cost_b_flow("network_simplex");

    //ABC_421_G("primal_dual");
    //ABC_421_G("cost_scaling");
    //ABC_421_G("network_simplex");

    //KUPC_2014_I("primal_dual");
    //KUPC_2014_I("cost_scaling");
    //KUPC_2014_I("network_simplex");

    //AOJ_2627("primal_dual");
    //AOJ_2627("cost_scaling");
    //AOJ_2627("network_simplex");

    //EducationalCodeforces80_F("primal_dual");
    //EducationalCodeforces80_F("cost_scaling");
    //EducationalCodeforces80_F("network_simplex");

    //UTPC2011_H("primal_dual");
    //UTPC2011_H("cost_scaling");
    //UTPC2011_H("network_simplex");

    //Codeforces826_G("primal_dual");
    //Codeforces826_G("cost_scaling");
    //Codeforces826_G("network_simplex");

    //ABC_393_G("primal_dual");
    //ABC_393_G("cost_scaling");
    ABC_393_G("network_simplex");
}