//
// 重み付き無向グラフの全域最小カット (by Stoer-Wagner 法, in O(N^3))
//
// reference:
//    https://dl.acm.org/doi/10.1145/263867.263872
//    https://cp-algorithms.com/graph/stoer_wagner_mincut.html
//


#include <bits/stdc++.h>
using namespace std;


// グラフを隣接行列形式で与える
template<class T> pair<T, vector<int>> stoer_wagner(vector<vector<T>> &G) {
    int V = (int)G.size();
    T best_cost = numeric_limits<T>::max()/2;
    vector<int> best_cut;
    vector<vector<int>> ids(V);
    for (int v = 0; v < V; v++) ids[v] = vector<int>(1, v);
    vector<T> w(V);
    vector<bool> exist(V, true), in_a(V);
    for (int iter = 0; iter < V - 1; iter++) {
        in_a.assign(V, false);
        w.assign(V, 0);
        for (int it = 0, prev; it < V - iter; it++) {
            int sel = -1;
            for (int v = 0; v < V; v++) {
                if (exist[v] && !in_a[v] && (sel == -1 || w[v] > w[sel])) {
                    sel = v;
                }
            }
            if (it == V - iter - 1) {
                if (w[sel] < best_cost) best_cost = w[sel], best_cut = ids[sel];
                ids[prev].insert(ids[prev].end(), ids[sel].begin(), ids[sel].end());
                for (int v = 0; v < V; v++) G[prev][v] = G[v][prev] += G[sel][v];
                exist[sel] = false;
            } else {
                in_a[sel] = true;
                for (int v = 0; v < V; v++) w[v] += G[sel][v];
                prev = sel;
            }
        }
    }
    sort(best_cut.begin(), best_cut.end());
    return {best_cost, best_cut};
}


//------------------------------//
// ナイーブ解と比較するための Dinic 法
//------------------------------//

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


//------------------------------//
// Examples
//------------------------------//

// ナイーブ解との比較
pair<int, vector<int>> verify
(int N, const vector<int> &U, const vector<int> &V, const vector<int> &W) {
    // Stoer-Wagner
    vector<vector<int>> G(N, vector<int>(N, 0));
    for (int i = 0; i < U.size(); i++) G[U[i]][V[i]] = G[V[i]][U[i]] = W[i];
    auto [cost, cut] = stoer_wagner(G);

    // Dinic
    int naive_res = 1LL << 29;
    FlowGraph<int> FG(N);
    for (int i = 0; i < U.size(); i++) FG.add_edge(U[i], V[i], W[i], W[i]);
    for (int s = 0; s < N; s++) for (int t = s+1; t < N; t++) {
        auto maxflow = Dinic(FG, s, t);
        FG.reset();
        naive_res = min(naive_res, maxflow);
    }
    assert(cost == naive_res);
    return {cost, cut};
}

unsigned int rand_int() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty; ty = tz; tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}
int rand_int(int minv, int maxv) {
    return rand_int() % (maxv - minv + 1) + minv;
}
void test_one_case(bool output = false) {
    int N = rand_int(2, 100);
    int M = rand_int(1, N * (N - 1) / 2);
    vector<int> U(M), V(M), W(M);
    set<pair<int,int>> already;
    for (int i = 0; i < M; i++) {
        U[i] = rand_int(0, N-2);
        V[i] = rand_int(U[i]+1, N-1);
        while (already.count({U[i], V[i]})) {
            U[i] = rand_int(0, N-2);
            V[i] = rand_int(U[i]+1, N-1);
        }
        already.insert({U[i], V[i]});
        W[i] = rand_int(1, 1000);
    }
    auto [cost, cut] = verify(N, U, V, W);
    if (output) {
        cout << "cost: " << cost << endl;
        cout << "cut: ";
        for (auto v : cut) cout << v << " ";
        cout << endl;
    }
}

void user_test() {
    for (int iter = 0; iter < 100; iter++) test_one_case(true);
}


int main() {
    user_test();
}