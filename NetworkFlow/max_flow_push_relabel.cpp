//
// max-flow (by Push-Relabel), in O(V^2√E)
//
// reference;
//   hitonanode: Maxflow (push-relabel, Goldberg & Tarjan) （Push-relabel による最大流）
//     https://hitonanode.github.io/cplib-cpp/flow/maxflow_pushrelabel.hpp
//
// verified:
//   典型アルゴリズム問題集 上級〜エキスパート編 E - 最大流
//     https://atcoder.jp/contests/pastbook2022/tasks/pastbook2022_e
//     https://atcoder.jp/contests/tessoku-book/tasks/tessoku_book_bp
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


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
// Examples
//------------------------------//

// 典型アルゴリズム問題集 上級〜エキスパート編 E - 最大流
void PAST_Max_Flow() {
    int V, E;
    cin >> V >> E;
    int s = 0, t = V - 1;
    FlowGraph<long long> G(V);
    for (int i = 0; i < E; ++i) {
        long long u, v, c;
        cin >> u >> v >> c, u--, v--;
        G.add_edge(u, v, c);
    }
    long long res = PushRelabel(G, s, t, false);
    cout << res << endl;
}


int main() {
    PAST_Max_Flow();
}