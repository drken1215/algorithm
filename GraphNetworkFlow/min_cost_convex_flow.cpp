//
// 最小凸費用流
//   各辺 e に、従来の「単位流量あたりのコスト」ではなく、凸関数 cost(f: 流量) を与える
//
// verified
//   AOJ Course GRL_6_B Network Flow - Minimum Cost Flow
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_6_B&lang=jp
//
//   JAG 模擬国内予選 2022 F - Nth Floor Dungeon and K Brave Men (AOJ 3297)
//     https://onlinejudge.u-aizu.ac.jp/challenges/sources/JAG/Prelim/3297?year=2022
//


#include <bits/stdc++.h>
using namespace std;


// edge class (for min convex-cost flow)
// cost-function must be convex
template<class FLOW, class COST> struct FlowConvexCostEdge {
    using FUNC = function<COST(FLOW)>;

    // core members
    int rev, from, to;
    FLOW cap, icap, flow;
    FUNC cost;
    
    // constructor
    FlowConvexCostEdge() {}
    FlowConvexCostEdge(int rev, int from, int to, FLOW cap, const FUNC &f)
        : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(0), cost(f) {
    }
    FlowConvexCostEdge(int rev, int from, int to, FLOW cap, FLOW rcap, const FUNC &f)
        : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(rcap), cost(f) {
    }
    constexpr void reset() { 
        flow -= icap - cap;
        cap = icap;
    }
    constexpr COST get_next_cost() const {
        return cost(flow + 1) - cost(flow);
    }
    friend constexpr COST get_next_cost(const FlowConvexCostEdge &e) {
        return e.get_next_cost();
    }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowConvexCostEdge& e) {
        return s << e.from << " -> " << e.to << " (" << e.cap << ", " << e.flow << ")";
    }
};

// graph class (for min-cost flow)
template<class FLOW, class COST> struct FlowConvexCostGraph {
    using FUNC = function<COST(FLOW)>;
    
    // core members
    vector<vector<FlowConvexCostEdge<FLOW, COST>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    vector<COST> pot; // pot[v] := potential (e.cost + pot[e.from] - pos[e.to] >= 0)
    bool include_negative_edge = false;
    
    // constructor
    FlowConvexCostGraph(int n = 0) : list(n), pot(n), include_negative_edge(false) { }
    void init(int n = 0) {
        list.clear(), list.resize(n);
        pos.clear();
        pot.assign(n, 0);
        include_negative_edge = false;
    }
    
    // getter
    vector<FlowConvexCostEdge<FLOW, COST>> &operator [] (int i) {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    const vector<FlowConvexCostEdge<FLOW, COST>> &operator [] (int i) const {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowConvexCostEdge<FLOW, COST> &get_rev_edge(const FlowConvexCostEdge<FLOW, COST> &e) {
        return list[e.to][e.rev];
    }
    FlowConvexCostEdge<FLOW, COST> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowConvexCostEdge<FLOW, COST> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowConvexCostEdge<FLOW, COST>> get_edges() const {
        vector<FlowConvexCostEdge<FLOW, COST>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowConvexCostEdge<FLOW, COST> &e : list[i]) e.reset();
        }
    }
    
    // add_edge
    void add_edge(int from, int to, FLOW cap, const FUNC &cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        auto negcost = [cost, cap](FLOW f) -> COST { return cost(cap - f); };
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowConvexCostEdge<FLOW, COST>(to_id, from, to, cap, 0, cost));
        list[to].push_back(FlowConvexCostEdge<FLOW, COST>(from_id, to, from, 0, cap, negcost));
        if (cost(1) - cost(0) < 0) include_negative_edge = true;
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
                if (pot[e.to] >= pot[cur] + e.get_next_cost()) pot[e.to] = pot[cur] + e.get_next_cost();
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
                if (pot[e.to] > pot[cur] + e.get_next_cost()) {
                    pot[e.to] = pot[cur] + e.get_next_cost();
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
    friend ostream& operator << (ostream& s, const FlowConvexCostGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// min-cost max-flow (<= limit_flow), slope ver.
template<class FLOW, class COST> vector<pair<FLOW, COST>>
MinCostFlowSlope(FlowConvexCostGraph<FLOW, COST> &G, int S, int T, FLOW limit_flow)
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
                if (e.cap >= FLOW(1)) {
                    COST add = e.get_next_cost() + G.pot[v] - G.pot[e.to];
                    if (dist[e.to] > dist[v] + add) {
                        dist[e.to] = dist[v] + add;
                        prevv[e.to] = v;
                        preve[e.to] = i;
                        que.emplace(dist[e.to], e.to);
                    }
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
        FLOW flow = 1;  // 1 ずつ流す
        COST cost = G.pot[T] - G.pot[S];
        for (int v = T; v != S; v = prevv[v]) {
            FlowConvexCostEdge<FLOW, COST> &e = G[prevv[v]][preve[v]];
            FlowConvexCostEdge<FLOW, COST> &re = G.get_rev_edge(e);
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
MinCostFlowSlope(FlowConvexCostGraph<FLOW, COST> &G, int S, int T)
{
    return MinCostFlowSlope(G, S, T, numeric_limits<FLOW>::max());
}

// min-cost max-flow (<= limit_flow)
template<class FLOW, class COST> pair<FLOW, COST>
MinCostFlow(FlowConvexCostGraph<FLOW, COST> &G, int S, int T, FLOW limit_flow)
{
    return MinCostFlowSlope(G, S, T, limit_flow).back();
}

// min-cost max-flow (<= limit_flow)
template<class FLOW, class COST> pair<FLOW, COST>
MinCostFlow(FlowConvexCostGraph<FLOW, COST> &G, int S, int T)
{
    return MinCostFlow(G, S, T, numeric_limits<FLOW>::max());
}

// Minimum Cost b-flow (come down to min-cost primal-dual, negative cycle is NG)
template<class FLOW, class COST> struct MinConvexCostBFlow {
    using FUNC = function<COST(FLOW)>;

    // Edge
    struct Edge {
        int from, to;
        FLOW lower_cap, upper_cap, flow;
        FUNC cost;

        // debug
        friend ostream& operator << (ostream& s, const Edge& e) {
            return s << e.from << "->" << e.to << '(' << e.lower_cap << '~' << e.upper_cap << ')';
        }
    };

    // inner values
    int V;
    vector<Edge> edges;
    vector<FLOW> lower_dss, upper_dss;  // demand (< 0) and supply (> 0)
    vector<COST> dual;
    FlowConvexCostGraph<FLOW, COST> G;

    // constructor
    MinConvexCostBFlow() {}
    MinConvexCostBFlow(int V) : V(V), lower_dss(V, 0), upper_dss(V, 0) {}

    // setter
    void add_edge(int from, int to, FLOW cap, const FUNC &cost) {
        assert(cap >= 0);
        edges.push_back({from, to, 0, cap, 0, cost});
    }
    void add_edge(int from, int to, FLOW lower_cap, FLOW upper_cap, const FUNC &cost) {
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
        G.init(V + 3);
        int s = V + 1, t = V + 2;
        for (const auto &e : edges) {
            dss[e.to] += e.lower_cap, dss[e.from] -= e.lower_cap;
            auto newcost = [e](FLOW f) -> COST { return e.cost(f + e.lower_cap); };
            G.add_edge(e.from, e.to, e.upper_cap - e.lower_cap, newcost);
        }

        // ds -> s, t
        FLOW ssum = 0, tsum = 0;
        auto zero = [](FLOW f) -> COST { return 0; };
        for (int i = 0; i < V + 1; i++) {
            if (dss[i] > 0) ssum += dss[i], G.add_edge(s, i, dss[i], zero);
            else if (dss[i] < 0) tsum -= dss[i], G.add_edge(i, t, -dss[i], zero);
        }

        // feasibility check
        if (ssum != tsum) return {false, COST(0)};
        
        // min-cost flow
        auto [maxflow, mincost] = MinCostFlow(G, s, t, ssum);
        if (maxflow < ssum) return {false, COST(0)};
        COST res = 0;
        for (int i = 0; i < (int)edges.size(); i++) {
            auto &e = edges[i];
            const auto &ge = G.get_edge(i);
            e.flow = e.upper_cap - ge.cap;
            res += e.cost(e.flow);
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
// Examples
//------------------------------//

// AOJ
void AOJ_Course_GRL_6_B() {
    int V, E, F;
    cin >> V >> E >> F;
    FlowConvexCostGraph<int, int> G(V);
    for (int i = 0; i < E; ++i) {
        int u, v, cap, cost;
        cin >> u >> v >> cap >> cost;
        auto func = [cost](int f) -> int { return cost * f; };
        G.add_edge(u, v, cap, func);
    }
    int s = 0, t = V-1;
    auto [max_flow, min_cost] = MinCostFlow(G, s, t, F);
    cout << (max_flow == F ? min_cost : -1) << endl;
}


// JAG 模擬国内予選 2022 F - Nth Floor Dungeon and K Brave Men (AOJ 3297)
void AOJ_3297() {
    long long N, M, K;
    while (cin >> N >> M >> K, N) {
        vector<long long> S(M), T(M), C(M), D(M);
        for (int i = 0; i < M; i++) cin >> S[i] >> T[i] >> C[i] >> D[i], S[i]--, T[i]--;

        MinConvexCostBFlow<long long, long long> G(N);
        G.set_ds(0, K), G.set_ds(N-1, -K);
        for (int i = 0; i < M; i++) {
            auto cost = [&C, &D, i](long long f) -> long long { 
                return C[i] * f + D[i] * (f - 1) * f / 2;
            };
            G.add_edge(S[i], T[i], 1, K, cost);
        }
        auto [flag, cost] = G.solve();
        cout << (flag ? cost : -1) << endl;
    }
}


int main() {
    //AOJ_Course_GRL_6_B();
    AOJ_3297();
}