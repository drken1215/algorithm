//
// b-flow 制約下で、s-t 間の最大フローを求める
//
// verified
//   AOJ 1615 - Gift Exchange Party (ICPC 国内予選 2016 H)
//     https://onlinejudge.u-aizu.ac.jp/problems/1615
//
//   AtCoder ABC 285 G - Tatami
//     https://atcoder.jp/contests/abc285/tasks/abc285_g
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
    FlowEdge() {}
    FlowEdge(int r, int f, int t, FLOW c) : rev(r), from(f), to(t), cap(c), icap(c), flow(0) {}
    void reset() { cap = icap, flow = 0; }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowEdge& E) {
        return s << E.from << "->" << E.to << '(' << E.flow << '/' << E.icap << ')';
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
    void change_edge(FlowEdge<FLOW> &e, FLOW new_cap, FLOW new_flow) {
        assert(0 <= new_flow && new_flow <= new_cap);
        FlowEdge<FLOW> &re = get_rev_edge(e);
        e.cap = new_cap - new_flow, e.icap = new_cap, e.flow = new_flow;
        re.cap = new_flow;
    }
    
    // add_edge
    void add_edge(int from, int to, FLOW cap, FLOW rcap = 0) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowEdge<FLOW>(to_id, from, to, cap));
        list[to].push_back(FlowEdge<FLOW>(from_id, to, from, rcap));
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

// Maximum b-flow
template<class FLOW> struct MaxBFlow {
    // Edge
    struct Edge {
        int from, to;
        FLOW lower_cap, upper_cap, flow;

        // debug
        friend ostream& operator << (ostream& s, const Edge& e) {
            return s << e.from << "->" << e.to 
            << '(' << e.lower_cap << '~' << e.upper_cap << ')';
        }
    };

    // inner values
    int V;
    vector<Edge> edges;
    vector<FLOW> lower_dss, upper_dss;  // demand (< 0) and supply (> 0)
    FlowGraph<FLOW> G;

    // constructor
    MaxBFlow() {}
    MaxBFlow(int V) : V(V), lower_dss(V, 0), upper_dss(V, 0) {}

    // setter
    void add_edge(int from, int to, FLOW cap) {
        assert(cap >= 0);
        edges.push_back({from, to, 0, cap, 0});
    }
    void add_edge(int from, int to, FLOW lower_cap, FLOW upper_cap) {
        assert(lower_cap <= upper_cap);
        edges.push_back({from, to, lower_cap, upper_cap, 0});
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
    pair<bool, FLOW> solve(int s, int t) {
        assert(0 <= s && s < V);
        assert(0 <= t && t < V);
        assert(s != t);

        // lower_ds, upper_ds -> strict ds
        int super = V;
        vector<FLOW> dss(V + 1, 0);
        for (int i = 0; i < V; i++) {
            if (lower_dss[i] == upper_dss[i]) dss[i] = lower_dss[i];
            else if (lower_dss[i] >= 0) {
                add_edge(super, i, lower_dss[i], upper_dss[i]);
            } else if (upper_dss[i] < 0) {
                add_edge(i, super, -upper_dss[i], -lower_dss[i]);
            } else {
                add_edge(super, i, upper_dss[i]);
                add_edge(i, super, -lower_dss[i]);
            }
        }

        // pre-flow lower_cap
        G.init(V + 3);
        for (const auto &e : edges) {
            dss[e.to] += e.lower_cap, dss[e.from] -= e.lower_cap;
            G.add_edge(e.from, e.to, e.upper_cap - e.lower_cap);
        }

        // ds -> s2, t2
        int s2 = V + 1, t2 = V + 2;
        FLOW ssum = 0, tsum = 0;
        for (int i = 0; i < V + 1; i++) {
            if (dss[i] > 0) ssum += dss[i], G.add_edge(s2, i, dss[i]);
            else if (dss[i] < 0) tsum -= dss[i], G.add_edge(i, t2, -dss[i]);
        }
        if (ssum != tsum) return {false, FLOW(-1)};

        // main solver
        FLOW a = Dinic(G, s2, t2);
        FLOW b = Dinic(G, s2, t);
        FLOW c = Dinic(G, s, t2); 
        FLOW d = Dinic(G, s, t);
        if (a + b != ssum || a + c != tsum) return {false, FLOW(-1)};
        else return {true, c + d};
    }
};


//------------------------------//
// Solver
//------------------------------//

// AOJ 1615 - Gift Exchange Party (ICPC 国内予選 2016 H)
void AOJ_1615() {
    int N, M;
    while (cin >> N >> M, N) {
        vector<int> U(M), V(M);
        for (int j = 0; j < M; j++) cin >> U[j] >> V[j], U[j]--, V[j]--;
        for (int gap = 0; gap <= N; gap++) {
            for (int mi = N - gap; mi >= 0; mi--) {
                int ma = mi + gap, s = N + M, t = s + 1;
                MaxBFlow<int> G(N + M + 2);
                for (int i = 0; i < N; i++) G.add_edge(s, i, mi, ma);
                for (int j = 0; j < M; j++) {
                    G.add_edge(j + N, t, 1);
                    G.add_edge(U[j], j + N, 1), G.add_edge(V[j], j + N, 1);
                }
                auto [flag, flow] = G.solve(s, t);
                if (flag && flow == M) {
                    cout << mi << " " << ma << endl;
                    goto END;
                }
            }
        }
        END:;
    }
}

// AtCoder ABC 285 G - Tatami
void ABC_285_G() {
    const vector<int> DX = {1, 0, -1, 0};
    const vector<int> DY = {0, 1, 0, -1};
    int H, W; cin >> H >> W;
    vector<string> C(H);
    for (int i = 0; i < H; i++) cin >> C[i];
    int V = H * W, s = V * 2, t = s + 1;
    MaxBFlow<int> G(V * 2 + 2);
    for (int i = 0; i < H; i++) for (int j = 0; j < W; j++) {
        int v = i * W + j;
        if (C[i][j] == '2') G.add_edge(s, v, 1, 1), G.add_edge(v+V, t, 1, 1);
        else G.add_edge(s, v, 0, 1), G.add_edge(v+V, t, 0, 1);
        for (int d = 0; d < 4; d++) {
            int i2 = i + DX[d], j2 = j + DY[d];
            if (i2 < 0 || i2 >= H || j2 < 0 || j2 >= W) continue;
            int v2 = i2 * W + j2;
            if (C[i][j] != '1' && C[i2][j2] != '1') G.add_edge(v, v2+V, 0, 1);
        }
    }
    auto [flag, flow] = G.solve(s, t);
    cout << (flag ? "Yes" : "No") << endl;
}


int main() {
    AOJ_1615();
    //ABC_285_G();
}