//
// min-cost flow (primal-dual)
//
// verified
//   典型アルゴリズム問題集 上級〜エキスパート編 F - 最小費用流
//     https://atcoder.jp/contests/pastbook2022/tasks/pastbook2022_f
//
//   AtCoder Library Practice Contest E - MinCostFlow
//     https://atcoder.jp/contests/practice2/tasks/practice2_e
//


#include <bits/stdc++.h>
using namespace std;


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
template<class FLOWTYPE, class COSTTYPE> vector<pair<FLOWTYPE, COSTTYPE>>
MinCostFlowSlope(FlowCostGraph<FLOWTYPE, COSTTYPE> &G, int S, int T, FLOWTYPE limit_flow)
{
    // result values
    FLOWTYPE cur_flow = 0;
    COSTTYPE cur_cost = 0, pre_cost = -1;
    vector<pair<FLOWTYPE, COSTTYPE>> res;
    res.emplace_back(cur_flow, cur_cost);
    
    // intermediate values
    vector<COSTTYPE> dual((int)G.size(), 0), dist((int)G.size());
    vector<int> prevv((int)G.size(), -1), preve((int)G.size(), -1);
    
    // dual
    auto dual_step = [&]() -> bool {
        priority_queue<pair<COSTTYPE,int>, vector<pair<COSTTYPE,int>>, greater<pair<COSTTYPE,int>>> que;
        que.push({0, S});
        dist.assign((int)G.size(), numeric_limits<COSTTYPE>::max());
        dist[S] = 0;
        while (!que.empty()) {
            auto [cur_cost, v] = que.top();
            que.pop();
            if (dist[v] < cur_cost) continue;
            for (int i = 0; i < (int)G[v].size(); ++i) {
                const auto &e = G[v][i];
                COSTTYPE new_cost = e.cost + dual[v] - dual[e.to];
                if (e.cap > 0 && dist[e.to] > dist[v] + new_cost) {
                    dist[e.to] = dist[v] + new_cost;
                    prevv[e.to] = v;
                    preve[e.to] = i;
                    que.push({dist[e.to], e.to});
                }
            }
        }
        if (dist[T] == numeric_limits<COSTTYPE>::max()) return false;
        for (int v = 0; v < (int)G.size(); ++v) {
            if (dist[T] == numeric_limits<COSTTYPE>::max()) continue;
            dual[v] -= dist[T] - dist[v];
        }
        return true;
    };
    
    // primal
    auto primal_step = [&]() -> void {
        FLOWTYPE flow = limit_flow - cur_flow;
        COSTTYPE cost = -dual[S];
        for (int v = T; v != S; v = prevv[v]) {
            flow = min(flow, G[prevv[v]][preve[v]].cap);
        }
        for (int v = T; v != S; v = prevv[v]) {
            FlowCostEdge<FLOWTYPE, COSTTYPE> &e = G[prevv[v]][preve[v]];
            FlowCostEdge<FLOWTYPE, COSTTYPE> &re = G.get_rev_edge(e);
            e.cap -= flow, e.flow += flow;
            re.cap += flow, re.flow -= flow;
        }
        cur_flow += flow;
        cur_cost += flow * cost;
        if (pre_cost == cost) res.pop_back();
        res.emplace_back(cur_flow, cur_cost);
        pre_cost = cur_cost;
    };
    
    // primal-dual
    while (cur_flow < limit_flow) {
        if (!dual_step()) break;
        primal_step();
    }
    return res;
}

// min-cost max-flow, slope ver.
template<class FLOWTYPE, class COSTTYPE> vector<pair<FLOWTYPE, COSTTYPE>>
MinCostFlowSlope(FlowCostGraph<FLOWTYPE, COSTTYPE> &G, int S, int T)
{
    return MinCostFlowSlope(G, S, T, numeric_limits<FLOWTYPE>::max());
}

// min-cost max-flow (<= limit_flow)
template<class FLOWTYPE, class COSTTYPE> pair<FLOWTYPE, COSTTYPE>
MinCostFlow(FlowCostGraph<FLOWTYPE, COSTTYPE> &G, int S, int T, FLOWTYPE limit_flow)
{
    return MinCostFlowSlope(G, S, T, limit_flow).back();
}

// min-cost max-flow (<= limit_flow)
template<class FLOWTYPE, class COSTTYPE> pair<FLOWTYPE, COSTTYPE>
MinCostFlow(FlowCostGraph<FLOWTYPE, COSTTYPE> &G, int S, int T)
{
    return MinCostFlow(G, S, T, numeric_limits<FLOWTYPE>::max());
}



//------------------------------//
// Examples
//------------------------------//

// 典型アルゴリズム問題集 上級〜エキスパート編 F - 最小費用流
void PAST_Min_Cost_Flow() {
    long long V, E, F;
    cin >> V >> E >> F;
    FlowCostGraph<long long, long long> G(V);
    for (int i = 0; i < E; ++i) {
        long long u, v, cap, cost;
        cin >> u >> v >> cap >> cost, u--, v--;
        G.add_edge(u, v, cap, cost);
    }
    long long s = 0, t = V-1;
    auto [max_flow, min_cost] = MinCostFlow(G, s, t, F);
    cout << (max_flow == F ? min_cost : -1) << endl;
}

// ACL practice E
void ACL_practice_E() {
    // 十分大きな値
    const long long B = 1LL<<40;
    
    // 入力
    int N, K;
    cin >> N >> K;
    vector<vector<long long>> A(N, vector<long long>(N));
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) cin >> A[i][j];
    
    // フローネットワークを作る
    // 行番号に対応する頂点を 0, 1, ..., N-1、列番号に対応する頂点を N, N+1, ..., 2N-1 とする
    // 超頂点の番号を S = 2N, T = 2N+1 とする
    FlowCostGraph<int, long long> G(N * 2 + 2);
    int S = N * 2, T = N * 2 + 1;
    
    // 行と列を結ぶ
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // 容量 1、コスト B - A[i][j]
            G.add_edge(i, j + N, 1, B - A[i][j]);
        }
    }
    
    // 超頂点
    for (int i = 0; i < N; ++i) {
        G.add_edge(S, i, K, 0);  // 容量 K, コスト 0
        G.add_edge(i + N, T, K, 0);  // 容量 K, コスト 0
    }
    
    // バイパス
    G.add_edge(S, T, N * K, B);
    
    // 流量 N * K の最小費用流を流す (最大流量も受け取るが N * K になることは分かっている)
    auto [max_flow, min_cost] = MinCostFlow(G, S, T, N * K);
    
    // 復元する
    vector<string> grid(N, string(N, '.'));
    const auto &edges = G.get_edges();
    for (const auto &e : edges) {
        // 超頂点が絡む辺や、フローの流れなかった辺はスキップ
        if (e.from == S || e.to == T || e.flow == 0) continue;
        
        // 行 e.from、列 e.to - N が選ばれる
        grid[e.from][e.to - N] = 'X';
    }
    
    // 出力
    cout << B * N * K - min_cost << endl;
    for (int i = 0; i < N; ++i) cout << grid[i] << endl;
}


int main() {
    PAST_Min_Cost_Flow();
    //ACL_practice_E();
}