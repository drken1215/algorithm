//
// min-cost flow (primal-dual, negative edges are ok but negative cycles are ng)
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


// edge class (for network-flow)
template<class FLOWTYPE, class COSTTYPE> struct FlowCostEdge {
    // core members
    int rev, from, to;
    FLOWTYPE cap, icap, flow;
    COSTTYPE cost;
    
    // constructor
    FlowCostEdge(int rev, int from, int to, FLOWTYPE cap, COSTTYPE cost)
    : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(0), cost(cost) {}
    void reset() { cap = icap, flow = 0; }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowCostEdge& e) {
        return s << e.from << " -> " << e.to
        << " (" << e.flow << "/" << e.icap << ", " << e.cost << ")";
    }
};

// graph class (for network-flow)
template<class FLOWTYPE, class COSTTYPE> struct FlowCostGraph {
    // core members
    vector<vector<FlowCostEdge<FLOWTYPE, COSTTYPE>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    
    // constructor
    FlowCostGraph(int n = 0) : list(n) { }
    void init(int n = 0) {
        list.assign(n, FlowCostEdge<FLOWTYPE, COSTTYPE>());
        pos.clear();
    }
    
    // getter
    vector<FlowCostEdge<FLOWTYPE, COSTTYPE>> &operator [] (int i) {
        return list[i];
    }
    const vector<FlowCostEdge<FLOWTYPE, COSTTYPE>> &operator [] (int i) const {
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowCostEdge<FLOWTYPE, COSTTYPE> &get_rev_edge
    (const FlowCostEdge<FLOWTYPE, COSTTYPE> &e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    FlowCostEdge<FLOWTYPE, COSTTYPE> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowCostEdge<FLOWTYPE, COSTTYPE> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowCostEdge<FLOWTYPE, COSTTYPE>> get_edges() const {
        vector<FlowCostEdge<FLOWTYPE, COSTTYPE>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowCostEdge<FLOWTYPE, COSTTYPE> &e : list[i]) e.reset();
        }
    }
    
    // add_edge
    void add_edge(int from, int to, FLOWTYPE cap, COSTTYPE cost) {
        pos.emplace_back(from, (int)list[from].size());
        list[from].push_back(FlowCostEdge<FLOWTYPE, COSTTYPE>
                             ((int)list[to].size(), from, to, cap, cost));
        list[to].push_back(FlowCostEdge<FLOWTYPE, COSTTYPE>
                           ((int)list[from].size() - 1, to, from, 0, -cost));
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
    vector<bool> seen((int)G.size(), false);
    vector<COSTTYPE> dist((int)G.size(), numeric_limits<COSTTYPE>::max());
    vector<int> prevv((int)G.size(), -1), preve((int)G.size(), -1);
    
    // dual
    auto dual_step = [&]() -> bool {
        seen.assign((int)G.size(), false);
        dist.assign((int)G.size(), numeric_limits<COSTTYPE>::max());
        seen[S] = true, dist[S] = 0;
        while (true) {
            bool update = false;
            for (int v = 0; v < (int)G.size(); ++v) {
                if (!seen[v]) continue;
                for (int i = 0; i < G[v].size(); ++i) {
                    const FlowCostEdge<FLOWTYPE, COSTTYPE> &e = G[v][i];
                    if (e.cap > 0 && (!seen[e.to] || dist[e.to] > dist[v] + e.cost)) {
                        dist[e.to] = dist[v] + e.cost;
                        prevv[e.to] = v;
                        preve[e.to] = i;
                        seen[e.to] = true;
                        update = true;
                    }
                }
            }
            if (!update) break;
        }
        return seen[T];
    };
    
    // primal
    auto primal_step = [&]() -> void {
        FLOWTYPE flow = limit_flow - cur_flow;
        COSTTYPE cost = dist[T];
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
    //AOJ_Course_GRL_6_B();
    ACL_practice_E();
}

