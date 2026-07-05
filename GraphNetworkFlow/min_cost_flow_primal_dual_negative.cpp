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
//   AtCoder ABC 214 H - Collecting (for DAG potential)
//     https://atcoder.jp/contests/abc214/tasks/abc214_h
//
//   AtCoder ABC 247 G - Dream Team (for SPFA potential and slope)
//     https://atcoder.jp/contests/abc247/tasks/abc247_g
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
        : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(0), cost(cost) {
    }
    FlowCostEdge(int rev, int from, int to, FLOW cap, FLOW rcap, COST cost)
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
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, 0, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, 0, cap, -cost));
        if (cost < 0) include_negative_edge = true;
    }
    void add_edge(int from, int to, FLOW cap, FLOW rcap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, rcap, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, rcap, cap, -cost));
        if (cost < 0) include_negative_edge = true;
    }
    void add_bidirected_edge(int from, int to, FLOW cap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        add_edge(from, to, cap, cap, cost);
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
        for (int v = 0; v < size(); v++) que.push(v), inque[v] = true;
        while (!que.empty()) {
            int cur = que.front();
            que.pop();
            inque[cur] = false;
            if (cnt[cur] > size()) return false;  // include negative-cycle
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

// AtCoder ABC 214 H - Collecting
// scc
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};
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
        if (record_reversed_edges) reversed_list[from].push_back(Edge(from, to, val));
        if (from != to) {
            if (record_edge_index) id[to][from] = (int)list[to].size(); 
            list[to].push_back(Edge(to, from, val));
            if (record_reversed_edges) reversed_list[to].push_back(Edge(to, from, val));
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
template<class T> struct SCC {
    // results
    vector<int> cmp;
    vector<vector<int>> groups;
    Graph<T> dag;
    
    // intermediate results
    vector<bool> seen;
    vector<int> vs, rvs;
    
    // constructor
    SCC() { }
    SCC(const Graph<T> &G) { 
        solve(G);
    }
    void init(const Graph<T> &G) { 
        solve(G);
    }

    // getter, compressed dag（v: node-id of compressed dag)
    int get_size(int v) const {
        return groups[v].size();
    }
    vector<int> get_group(int v) const {
        return groups[v];
    }

    // solver
    void dfs(const Graph<T> &G, int v) {
        seen[v] = true;
        for (const auto &e : G[v]) if (!seen[e.to]) dfs(G, e.to);
        vs.push_back(v);
    }
    void rdfs(const Graph<T> &G, int v, int k) {
        seen[v] = true;
        cmp[v] = k;
        for (const auto &e : G.get_rev_edges(v)) if (!seen[e.to]) rdfs(G, e.to, k);
        rvs.push_back(v);
    }
    void reconstruct(const Graph<T> &G) {
        dag.init((int)groups.size());
        set<pair<int,int>> new_edges;
        for (int i = 0; i < (int)G.size(); ++i) {
            int u = cmp[i];
            for (const auto &e : G[i]) {
                int v = cmp[e.to];
                if (u == v) continue;
                if (!new_edges.count({u, v})) {
                    dag.add_edge(u, v);
                    new_edges.insert({u, v});
                }
            }
        }
    }
    void solve(const Graph<T> &G) {
        // first dfs
        seen.assign((int)G.size(), false);
        vs.clear();
        for (int v = 0; v < (int)G.size(); ++v) if (!seen[v]) dfs(G, v);

        // back dfs
        int k = 0;
        groups.clear();
        seen.assign((int)G.size(), false);
        cmp.assign((int)G.size(), -1);
        for (int i = (int)G.size()-1; i >= 0; --i) {
            if (!seen[vs[i]]) {
                rvs.clear();
                rdfs(G, vs[i], k++);
                groups.push_back(rvs);
            }
        }
        reconstruct(G);
    }
};
void ABC_214_H() {
    long long N, M, K;
    cin >> N >> M >> K;
    Graph G(N, true);
    for (int i = 0; i < M; i++) {
        int A, B;
        cin >> A >> B, A--, B--;
        G.add_edge(A, B);
    }
    vector<long long> X(N);
    for (int v = 0; v < N; v++) cin >> X[v];
    SCC scc(G);
    auto cmp = scc.cmp;
    auto dag = scc.dag;
    long long V = dag.size();
    vector<long long> W(V, 0);
    for (int v = 0; v < N; v++) W[cmp[v]] += X[v];

    const long long INF = 1LL<<50;
    FlowCostGraph<long long, long long> FG(V*2+1);
    long long t = V*2;
    for (int v = 0; v < V; v++) {
        FG.add_edge(v, v+V, 1, -W[v]);
        FG.add_edge(v, v+V, K-1, 0);
        for (auto e : dag[v]) FG.add_edge(e.from+V, e.to, K, 0);
        FG.add_edge(v+V, t, K, 0);
    }
    auto [flow, cost] = MinCostFlow(FG, cmp[0], t, K);
    cout << -cost << endl;
}

// AtCoder ABC 247 G - Dream Team (for SPFA potential)
void ABC_247_G() {
    long long N, M = 200, INF = 1LL<<29;

    cin >> N;
    vector<long long> A(N), B(N), C(N);
    for (int i = 0; i < N; i++) cin >> A[i] >> B[i] >> C[i], A[i]--, B[i]--;

    FlowCostGraph<long long, long long> FG(M * 2 + 2);
    long long s = M * 2, t = s + 1;
    for (int i = 0; i < M; i++) FG.add_edge(s, i, 1, 0), FG.add_edge(i+M, t, 1, 0);
    for (int i = 0; i < N; i++) FG.add_edge(A[i], B[i]+M, 1, -C[i]);

    auto slope = MinCostFlowSlope(FG, s, t);
    long long K = slope.back().first;
    vector<long long> res(K+1, 0);
    for (int i = 0; i < (int)slope.size()-1; i++) {
        auto [x1, y1] = slope[i];
        auto [x2, y2] = slope[i+1];
        for (long long x = x1; x <= x2; x++) {
            res[x] = (y2 - y1) / (x2 - x1) * (x - x1) + y1;
        }
    }
    cout << K << endl;
    for (int k = 1; k <= K; k++) cout << -res[k] << endl;
}


int main() {
    //AOJ_Course_GRL_6_B();
    //ACL_practice_E();
    //ABC_214_H();
    ABC_247_G();
}