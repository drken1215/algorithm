//
// Gomory-Hu 木 (in O(N f(N, M)), f(N, M): 最大フローを求める計算量）
//   無向グラフのすべての頂点対の最小カットの値が、木上の対応する頂点間のパス上にある辺重みの最小値に一致するような木
//
// Example
//   Codeforces Round 200 (Div. 1) E. Pumping Stations
//     https://codeforces.com/contest/343/problem/E
//
//   パ研合宿2024　第1日「SpeedRun」 R - Maximum Water Flow
//     https://atcoder.jp/contests/pakencamp-2024-day1/tasks/pakencamp_2024_day1_r
//


#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Max Flow (subroutine)
//------------------------------//

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

// Push-Relabel
// we can skip 2nd phase if we should know only about maxflow and residual graph
template<class FLOW> FLOW PushRelabel
(FlowGraph<FLOW> &G, int s, int t, FLOW limit_flow, bool do_2nd_phase = false) {
    assert(0 <= s && s < (int)G.size());
    assert(0 <= t && t < (int)G.size());
    assert(s != t);
    const int GlobalRelabelRreq = 5;
    const bool UseGapRelabeling = true;
    struct PushQueue {
        vector<pair<int, int>> even, odd;
        int num_even, num_odd;
        void init(int N) { even.resize(N), odd.resize(N), num_even = num_odd = 0; }
        void clear() { num_even = num_odd = 0; }
        int size() const { return num_even + num_odd; }
        bool empty() const { return size() == 0; }
        int highest() const {
            int a = (num_even > 0 ? even[num_even - 1].second : -1);
            int b = (num_odd > 0 ? odd[num_odd - 1].second : -1);
            return (a > b ? a : b);
        }
        void push(int v, int h) {
            if (h & 1) odd[num_odd++] = {v, h};
            else even[num_even++] = {v, h};
        }
        int pop() {
            if (num_even == 0 || (num_odd > 0 && odd[num_odd - 1].second > even[num_even - 1].second)) {
                return odd[--num_odd].first;
            } else {
                return even[--num_even].first;
            }
        }
    } push_que;

    int gap, N = (int)G.size();
    vector<int> dist, dcnt;
    vector<FLOW> excess;

    // heuristics
    auto global_relabeling = [&](int t) -> void {
        push_que.clear();
        if (UseGapRelabeling) gap = 1, dcnt.assign(N + 1, 0);
        dist.assign(N, N);
        dist[t] = 0;
        static vector<int> que;
        if (que.empty()) que.resize(N);
        que[0] = t;
        int qb = 0, qe = 1;
        while (qb < qe) {
            int now = que[qb++];
            if (UseGapRelabeling) gap = dist[now] + 1, dcnt[dist[now]]++;
            if (excess[now] > 0) push_que.push(now, dist[now]);
            for (const auto &e : G[now]) {
                if (G.get_rev_edge(e).cap > 0 && dist[e.to] == N) {
                    dist[e.to] = dist[now] + 1;
                    while ((int)que.size() <= qe) que.emplace_back(0);
                    que[qe++] = e.to;
                }
            }
        }
    };

    // push
    auto push = [&](int v, FlowEdge<FLOW> &e) -> void {
        auto &re = G.get_rev_edge(e);
        FLOW delta = e.cap < excess[v] ? e.cap : excess[v];
        excess[v] -= delta, e.cap -= delta, e.flow += delta;
        excess[e.to] += delta, re.cap += delta, re.flow -= delta;
        if (excess[e.to] > 0 && excess[e.to] <= delta) {
            if (!UseGapRelabeling || dist[e.to] <= gap) push_que.push(e.to, dist[e.to]);
        }
    };

    // run
    auto run = [&](int t) -> void {
        global_relabeling(t);
        int tick = (int)G.pos.size() * GlobalRelabelRreq;
        while (!push_que.empty()) {
            int v = push_que.pop();
            if (UseGapRelabeling && dist[v] > gap) continue;
            int dnex = N * 2 - 1;
            for (auto &e : G[v]) {
                if (e.cap <= 0) continue;
                if (dist[e.to] == dist[v] - 1) {
                    push(v, e);
                    if (excess[v] <= 0) break;
                } else {
                    if (dist[e.to] + 1 < dnex) dnex = dist[e.to] + 1;
                }
            }
            if (excess[v] > 0) {
                if (UseGapRelabeling) {
                    if (dnex != dist[v] && dcnt[dist[v]] == 1 && dist[v] < gap) gap = dist[v];
                    if (dnex == gap) gap++;
                    while (push_que.highest() > gap) push_que.pop();
                    if (dnex > gap) dnex = N;
                    if (dist[v] != dnex) dcnt[dist[v]]--, dcnt[dnex]++;
                }
                dist[v] = dnex;
                if (!UseGapRelabeling || dist[v] < gap) push_que.push(v, dist[v]);
            }
            if (GlobalRelabelRreq && --tick == 0) {
                tick = (int)G.pos.size() * GlobalRelabelRreq;
                global_relabeling(t);
            }
        }
    };

    // 1st phase: find preflow
    excess.assign(N, 0), dist.assign(N, 0);
    excess[s] += limit_flow, excess[t] -= limit_flow;
    dist[s] = N;
    if (UseGapRelabeling) gap = 1, dcnt.assign(N + 1, 0), dcnt[0] = N - 1;
    push_que.init(N);
    for (auto &e : G[s]) push(s, e);
    run(t);
    FLOW res = excess[t] + limit_flow;

    // 2nd phase: convert preflow into flow
    if (do_2nd_phase) {
        excess[s] += excess[t], excess[t] = 0;
        global_relabeling(s);
        run(s);
        assert(excess == vector<FLOW>(N, 0));
    }
    return res;
}

template<class FLOW> FLOW PushRelabel
(FlowGraph<FLOW> &G, int s, int t, bool do_2nd_phase = false) {
    return PushRelabel(G, s, t, numeric_limits<FLOW>::max(), do_2nd_phase);
}


//------------------------------//
// Gomory-Hu 木
//------------------------------//å

// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
    Edge(const Edge&) = default;
    Edge& operator = (const Edge&) = default;
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
template<class T = long long> struct Graph {
    int V;
    bool record_reversed_edges = false, record_edge_index = false;
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    vector<unordered_map<int, int>> id;  // id[v][w] := the index of node w in G[v]

    // constructors
    Graph(int n = 0, bool rre = false, bool rei = false) {
        init(n, rre, rei);
    }
    void init(int n = 0, bool rre = false, bool rei = false) {
        V = n, record_reversed_edges = rre, record_edge_index = rei;
        list.assign(n, vector<Edge<T>>());
        if (record_reversed_edges) reversed_list.assign(n, vector<Edge<T>>());
        if (record_edge_index) id.assign(n, unordered_map<int, int>());
    }
    Graph(const Graph&) = default;
    Graph& operator = (const Graph&) = default;

    // getters
    vector<Edge<T>> &operator [] (int i) { return list[i]; }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    constexpr size_t size() const { return list.size(); }
    constexpr void clear() { V = 0; list.clear(); }
    constexpr void resize(int n) { V = n; list.resize(n); }
    const vector<Edge<T>> &get_rev_edges(int i) const { 
        assert(record_reversed_edges);
        return reversed_list[i];
    }
    Edge<T> &get_edge(int u, int v) {
        assert(record_edge_index);
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u][v] >= 0 && id[u][v] < list[u].size());
        return list[u][id[u][v]];
    }
    const Edge<T> &get_edge(int u, int v) const {
        assert(record_edge_index);
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u].at(v) >= 0 && id[u].at(v) < list[u].size());
        return list[u][id[u].at(v)];
    }

    // add edge
    void add_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size(); 
        list[from].push_back(Edge(from, to, val));
        if (record_reversed_edges) reversed_list[to].push_back(Edge(to, from, val));
    }
    void add_bidirected_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size();
        list[from].push_back(Edge(from, to, val));
        if (record_reversed_edges) reversed_list[to].push_back(Edge(to, from, val));
        if (from != to) {
            if (record_edge_index) id[to][from] = (int)list[to].size(); 
            list[to].push_back(Edge(to, from, val));
            if (record_reversed_edges) reversed_list[from].push_back(Edge(from, to, val));
        }
    }

    // input (only tree-case)
    friend istream& operator >> (istream &is, Graph &G) {
        for (int i = 0; i < G.V - 1; i++) {
            int u, v;
            is >> u >> v, u--, v--;
            G.add_bidirected_edge(u, v);
        }
        return is;
    }

    // output
    friend ostream &operator << (ostream &os, const Graph &G) {
        os << endl;
        for (int i = 0; i < (int)G.size(); ++i) {
            os << i << " -> ";
            for (int j = 0; j < (int)G[i].size(); j++) {
                if (j) os << ", ";
                os << G[i][j].to << "(" << G[i][j].val << ")";
            }
            os << endl;
        }
        return os;
    }
};

// Gomory-Hu 木
template<class FLOW> Graph<FLOW> GomoryHuTree(const FlowGraph<FLOW> &G) {
    const FLOW INF = numeric_limits<FLOW>::max() / 2;
    vector<FlowEdge<FLOW>> res;
    int N = (int)G.size(), m = 1;
    const auto &es = G.get_edges();
    vector<vector<int>> vs(N);
    for (int i = 0; i < N; i++) vs[0].push_back(i);
    for (int i = 0; i < N; i++) {
        while (vs[i].size() > 1) {
            FlowGraph<FLOW> SG(N);
            for (int j = 0; j < i; j++) for (int k = 0; k < (int)vs[j].size() - 1; k++) {
                SG.add_bidirected_edge(vs[j][k], vs[j][k+1], INF);
            }
            for (const auto &e : es) SG.add_bidirected_edge(e.from, e.to, e.cap);
            for (const auto &e : res) if (e.from != i && e.to != i) {
                SG.add_bidirected_edge(vs[e.from][0], vs[e.to][0], INF);
            }
            int s = vs[i][0], t = vs[i][1];
            int flow = Dinic(SG, s, t);
            vector<int> cut = SG.find_cut(s, t);
            vector<int> vs1, vs2;
            for (auto v : vs[i]) {
                if (cut[v] == 1) vs1.push_back(v);
                else vs2.push_back(v);
            }
            for (auto &e : res) {
                if (e.to == i) swap(e.from, e.to);
                if (e.from == i && cut[vs[e.to][0]] != 1) e.from = m;
            }
            res.push_back(FlowEdge<FLOW>(-1, i, m, flow, flow));
            vs[i] = vs1;
            vs[m++] = vs2;
        }
    }
    Graph<FLOW> resG(N);
    for (auto &e : res) {
        e.from = vs[e.from][0], e.to = vs[e.to][0];
        resG.add_bidirected_edge(e.from, e.to, e.cap);
    }
    return resG;
}


//------------------------------//
// Examples
//------------------------------//

// Codeforces Round 200 (Div. 1) E. Pumping Stations
void Codeforces200_E() {
    const int INF = 1 << 29;
    int N, M, u, v, c;
    cin >> N >> M;
    FlowGraph<int> FG(N);
    for (int i = 0; i < M; i++) {
        cin >> u >> v >> c, u--, v--;
        FG.add_bidirected_edge(u, v, c);
    }
    auto G = GomoryHuTree(FG);

    int res = 0, ma = -1, start = 0;
    for (int v = 0; v < N; v++) for (auto e : G[v]) {
        res += e.val;
        if (ma < e.val) ma = e.val, start = v;
    }
    res /= 2;
    vector<int> path;
    priority_queue<pair<int,int>> que;
    que.push({INF, start});
    vector<bool> seen(N, false);
    while (!que.empty()) {
        auto [cur, v] = que.top();
        que.pop();
        if (seen[v]) continue;
        seen[v] = true;
        path.push_back(v);
        for (auto e : G[v]) que.push({e.val, e.to});
    }
    cout << res << endl;
    for (auto v : path) cout << v+1 << " ";
    cout << endl;
}


// パ研合宿2024　第1日「SpeedRun」 R - Maximum Water Flow
template<class VAL> struct Hungarian {
    // inner values
    int L, R;  // left size, right size
    vector<vector<VAL>> G;

    // results
    vector<int> lr, rl;
    vector<VAL> dualL, dualR;
    
    // constructor
    explicit Hungarian(const vector<vector<VAL>> &G_) 
        : L((int)G_.size()), R((int)G_[0].size()), G(G_) { 
        assert(L <= R);
    }

    // getter
    vector<int> get_lr() const { return lr; };
    vector<int> get_rl() const { return rl; };
    pair<vector<VAL>, vector<VAL>> get_dual() const { return {dualL, dualR}; };
    
    // solver
    VAL solve() {
        lr.assign(L, -1), rl.assign(R, -1);
        dualL.assign(L, VAL(0)), dualR.assign(R, VAL(0));
        vector<VAL> dist(R);
        vector<int> index(R), prev(R);
        iota(index.begin(), index.end(), 0);

        auto calc_residue = [&](int i, int j) -> VAL { return G[i][j] - dualR[j]; };
        for (int f = 0; f < L; f++) {
            for (int j = 0; j < R; j++) dist[j] = calc_residue(f, j), prev[j] = f;
            VAL w = 0;
            int j = 0, l = 0;
            for (int s = 0, t = 0;;) {
                if (s == t) {
                    l = s, w = dist[index[t++]];
                    for (int k = t; k < R; k++) {
                        j = index[k];
                        if (dist[j] <= w) {
                            if (dist[j] < w) t = s, w = dist[j];
                            index[k] = index[t], index[t++] = j;
                        }
                    }
                    for (int k = s; k < t; k++) {
                        j = index[k];
                        if (rl[j] < 0) goto augment;
                    }
                }
                int q = index[s++], i = rl[q];
                for (int k = t; k < R; k++) {
                    j = index[k];
                    VAL h = calc_residue(i, j) - calc_residue(i, q) + w;
                    if (h < dist[j]) {
                        dist[j] = h, prev[j] = i;
                        if (h == w) {
                            if (rl[j] < 0) goto augment;
                            index[k] = index[t], index[t++] = j;
                        }
                    }
                }
            }
        augment:
            for (int k = 0; k < l; k++) dualR[index[k]] += dist[index[k]] - w;
            int i = 0;
            do {
                rl[j] = i = prev[j];
                swap(j, lr[i]);
            } while (i != f);
        }
        VAL res = 0;
        for (int i = 0; i < L; i++) {
            res += G[i][lr[i]];
            dualL[i] = G[i][lr[i]] - dualR[lr[i]];
        }
        return res;
    }
};
void Paken2024Day1_R() {
    const long long INF = 1LL << 50;
    long long N, M, u, v, w;
    cin >> N >> M;
    FlowGraph<long long> FG(N);
    for (int i = 0; i < M; i++) {
        cin >> u >> v >> w, u--, v--;
        FG.add_bidirected_edge(u, v, w);
    }
    auto G = GomoryHuTree(FG);
    auto calc = [&](auto &&calc, int v, int t, int p = -1) -> long long {
        if (v == t) return INF;
        for (auto e : G[v]) {
            if (e.to == p) continue;
            long long x = calc(calc, e.to, t, v);
            if (x != -1) return min(e.val, x);
        }
        return -1;
    };
    vector<vector<long long>> cost(N, vector<long long>(N));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        if (i == j) cost[i][j] = INF;
        else cost[i][j] = -calc(calc, i, j);
    }
    Hungarian<long long> opt(cost);
    cout << -opt.solve() << endl;
}


int main() {
    //Codeforces200_E(); 
    Paken2024Day1_R();
}