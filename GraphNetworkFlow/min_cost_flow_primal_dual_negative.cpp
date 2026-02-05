//
// min-cost flow (primal-dual, negative edges are ok but negative cycles are ng)
//   負辺がある場合には、ポテンシャル法で解消している
//     1: DAG であるとき、DAG 上の DP でポテンシャルを求める
//     2: DAG ではないが負閉路を含まないとき、SPFA でポテンシャルを求める
//
// verified
//   AOJ Course GRL_6_B Network Flow - Minimum Cost Flow
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_6_B&lang=jp
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
    FlowCostEdge() {}
    FlowCostEdge(int rev, int from, int to, FLOW cap, COST cost)
    : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(0), cost(cost) {}
    void reset() { cap = icap, flow = 0; }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowCostEdge& e) {
        return s << e.from << " -> " << e.to << " (" << e.flow << "/" << e.icap << ", " << e.cost << ")";
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
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, 0, -cost));
        if (cost < 0) include_negative_edge = true;
    }
    void add_edge(int from, int to, FLOW cap, FLOW rcap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, rcap, -cost));
        if (cost < 0) include_negative_edge = true;
    }

    // find initial potential (to resolve initial negative-edge)
    // pot[v] := potential (e.cost + pot[e.from] - pos[e.to] >= 0)
    bool init_potential_dag() {
        vector<int> deg(size(), 0), st;
        for (int v = 0; v < size(); v++) for (const auto &e : list[v]) deg[e.to] += (e.cap > 0);
        st.reserve(size());
        for (int v = 0; v < size(); v++) if (!deg[v]) st.emplace_back(v);
        for (int i = 0; i < size(); i++) {
            if (st.size() == i) return false;  // not DAG
            int cur = st[i];
            for (const auto &e : list[cur]) {
                if (!e.cap) continue;
                deg[e.to]--;
                if (deg[e.to] == 0) st.emplace_back(e.to);
                if (pot[e.to] >= pot[cur] + e.cost) pot[e.to] = pot[cur] + e.cost;
            }
        }
        return true;
    }
    bool init_potential_spfa() {
        queue<int> que;
        vector<bool> inque(size(), true);
        vector<int> cnt(size(), 0);
        for (int v = 0; v < size(); v++) que.push(v);
        while (!que.empty()) {
            int cur = que.front();
            que.pop();
            inque[cur] = false;
            if (cnt[cur] > size()) return false;  // include negative-cycle
            cnt[cur]++;
            for (const auto &e : list[cur]) {
                if (!e.cap) continue;
                if (pot[e.to] > pot[cur] + e.cost) {
                    pot[e.to] = pot[cur] + e.cost;
                    if (!inque[e.to]) inque[e.to] = true, que.push(e.to);
                }
            }
        }
        return true;
    }
    bool init_potential() {
        if (!include_negative_edge) return true;
        return init_potential_dag() || init_potential_spfa();
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
    COST cur_cost = 0, pre_cost = -1;
    vector<pair<FLOW, COST>> res;
    res.emplace_back(cur_flow, cur_cost);
    
    // intermediate values
    vector<COST> dist((int)G.size(), numeric_limits<COST>::max());
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
                if (e.cap && dist[e.to] > dist[v] + add) {
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
        pre_cost = cur_cost;
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



//------------------------------//
// Examples
//------------------------------//

// AOJ
void AOJ_Course_GRL_6_B() {
    int V, E, F;
    cin >> V >> E >> F;
    FlowCostGraph<int, int> G(V);
    for (int i = 0; i < E; ++i) {
        int u, v, cap, cost;
        cin >> u >> v >> cap >> cost;
        G.add_edge(u, v, cap, cost);
    }
    int s = 0, t = V-1;
    auto [max_flow, min_cost] = MinCostFlow(G, s, t, F);
    cout << (max_flow == F ? min_cost : -1) << endl;
}

// ACL practice E
void ACL_practice_E() {
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
            // 容量 1、コスト -A[i][j]
            G.add_edge(i, j + N, 1, -A[i][j]);
        }
    }
    
    // 超頂点
    for (int i = 0; i < N; ++i) {
        G.add_edge(S, i, K, 0);  // 容量 K, コスト 0
        G.add_edge(i + N, T, K, 0);  // 容量 K, コスト 0
    }
    
    // バイパス
    G.add_edge(S, T, N * K, 0);
    
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
    cout << -min_cost << endl;
    for (int i = 0; i < N; ++i) cout << grid[i] << endl;
}


int main() {
    AOJ_Course_GRL_6_B();
    //ACL_practice_E();
}

