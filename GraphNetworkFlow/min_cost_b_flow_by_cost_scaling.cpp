//
// 最小費用 b-flow：最小費用循環流に帰着して、Cost-Scaling で解く
//   max-flow 計算を O(f(V, E)) として、O(f(V, E)log(VC))
//
// verified
//   Yosupo Library Checker - Minimum Cost b-flow
//     https://judge.yosupo.jp/problem/min_cost_b_flow
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Utility
//------------------------------//

template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }

using pint = pair<int, int>;
using pll = pair<long long, long long>;
using tint = array<int, 3>;
using tll = array<long long, 3>;
using fint = array<int, 4>;
using fll = array<long long, 4>;
using qint = array<int, 5>;
using qll = array<long long, 5>;
using vint = vector<int>;
using vll = vector<long long>;
using ll = long long;
using u32 = unsigned int;
using u64 = unsigned long long;
using i128 = __int128_t;
using u128 = __uint128_t;
template <class T>
using min_priority_queue = priority_queue<T, vector<T>, greater<T>>;

#define REP(i, a) for (long long i = 0; i < (long long)(a); i++)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); i++)
#define RREP(i, a) for (long long i = (a)-1; i >= (long long)(0); --i)
#define RREP2(i, a, b) for (long long i = (b)-1; i >= (long long)(a); --i)
#define EB emplace_back
#define PB push_back
#define MP make_pair
#define MT make_tuple
#define FI first
#define SE second
#define ALL(x) x.begin(), x.end()
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl

// debug stream
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, array<T, 3> P)
{ return s << '<' << P[0] << ", " << P[1] << ", " << P[2] << '>'; }
template<class T> ostream& operator << (ostream &s, array<T, 4> P)
{ return s << '<' << P[0] << ", " << P[1] << ", " << P[2] << ", " << P[3] << '>'; }
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, deque<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, vector<vector<T> > P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, set<T> P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, multiset<T> P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, unordered_set<T> P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T1, class T2> ostream& operator << (ostream &s, map<T1,T2> P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
template<class T1, class T2> ostream& operator << (ostream &s, unordered_map<T1,T2> P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }

// int 128
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


//------------------------------//
// Flow
//------------------------------//

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
            if (level[v] >= level[e.to] || e.cap == 0) continue;
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
            if (!flow) break;
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
    
    // constructor
    FlowCostGraph(int n = 0) : list(n) { }
    void init(int n = 0) {
        list.clear(), list.resize(n);
        pos.clear();
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
    }
    void add_edge(int from, int to, FLOW cap, FLOW rcap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, rcap, -cost));
    }

    // debug
    friend ostream& operator << (ostream& s, const FlowCostGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// Min Cost Circulation Flow by Cost-Scaling 
template<class FLOW, class COST> COST MinCostCirculation(FlowCostGraph<FLOW, COST> &G) {
    COST eps = 0;
    vector<FLOW> balance(G.size(), 0);
    vector<COST> price(G.size(), 0);
    
    auto newcost = [&](const FlowCostEdge<FLOW, COST> &e) -> COST {
        return e.cost * (COST)G.size() - price[e.from] + price[e.to];
    };

    auto ConstructGaux = [&]() -> void {
        vector<bool> visited(G.size(), false);
        auto dfs = [&](auto &&dfs, int v) -> void {
            visited[v] = true;
            for (int i = 0; i < G[v].size(); ++i) {
                FlowCostEdge<FLOW, COST> &e = G[v][i];
                if (e.cap > 0 && !visited[e.to] && newcost(e) < 0) dfs(dfs, e.to);
            }
        };
        for (int v = 0; v < G.size(); ++v) if (balance[v] > 0) dfs(dfs, v);
        for (int v = 0; v < G.size(); ++v) if (visited[v]) price[v] += eps;
    };

    auto augment_blocking_flow = [&]() -> bool {
        vector<COST> iter(G.size(), 0);
        auto augment = [&](auto &&augment, int v, FLOW flow) -> FLOW {
            if (balance[v] < 0) {
                FLOW dif = min(flow, -balance[v]);
                balance[v] += dif;
                return dif;
            }
            for (; iter[v] < G[v].size(); iter[v]++) {
                auto &e = G[v][iter[v]], &re = G.get_rev_edge(e);
                if (e.cap > 0 && newcost(e) < 0) {
                    FLOW dif = augment(augment, e.to, min(flow, e.cap));
                    if (dif > 0) {
                        e.cap -= dif, e.flow += dif;
                        re.cap += dif, re.flow -= dif;
                        return dif;
                    }
                }
            }
            return 0;
        };
        bool finish = true;
        for (int v = 0; v < G.size(); ++v) {
            FLOW flow;
            while (balance[v] > 0 && (flow = augment(augment, v, balance[v])) > 0)
                balance[v] -= flow;
            if (balance[v] > 0) finish = false;
        }
        if (finish) return true;
        else return false;
    };

    for (int v = 0; v < G.size(); ++v) {
        for (int i = 0; i < G[v].size(); ++i) {
            FlowCostEdge<FLOW, COST> &e = G[v][i];
            if (e.cap > 0) eps = max(eps, -e.cost * (COST)G.size());
        }
        price[v] = 0;
    }
    while (eps > 1) {
        eps /= 2;
        for (int v = 0; v < G.size(); ++v) {
            for (int i = 0; i < G[v].size(); ++i) {
                auto &e = G[v][i], &re = G.get_rev_edge(e);
                if (e.cap > 0 && newcost(e) < 0) {
                    FLOW flow = e.cap;
                    balance[e.from] -= flow, balance[e.to] += flow;
                    e.cap -= flow, e.flow += flow;
                    re.cap += flow, re.flow -= flow;
                }
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

// Minimum Cost b-flow (come down to min-cost circulation)
template<class FLOW, class COST> struct MinCostBFlow {
    // Edge
    struct Edge {
        int from, to;
        FLOW lower_cap, upper_cap, flow;
        COST cost;
    };

    // inner values
    int V;
    vector<Edge> edges;
    vector<FLOW> dss;  // demand (< 0) and supply (> 0)
    vector<COST> dual;
    FlowCostGraph<FLOW, COST> G;

    // constructor
    MinCostBFlow() {}
    MinCostBFlow(int V) : V(V), dss(V) {}

    // setter
    void add_edge(int from, int to, FLOW lower_cap, FLOW upper_cap, COST cost) {
        assert(lower_cap <= upper_cap);
        edges.push_back({from, to, lower_cap, upper_cap, 0, cost});
    }
    void add_ds(int v, FLOW ds) {
        assert(0 <= v && v < V);
        dss[v] += ds;
    }

    // solver
    pair<bool, COST> solve() {
        // feasibility check
        FlowGraph<FLOW> sg(V + 2);
        int s = V, t = s + 1;
        for (const auto &e : edges) {
            dss[e.to] += e.lower_cap, dss[e.from] -= e.lower_cap;
            sg.add_edge(e.from, e.to, e.upper_cap - e.lower_cap);
        }
        FLOW ssum = 0, tsum = 0;
        for (int i = 0; i < V; i++) {
            if (dss[i] > 0) ssum += dss[i], sg.add_edge(s, i, dss[i]);
            else if (dss[i] < 0) tsum -= dss[i], sg.add_edge(i, t, -dss[i]);
        }
        if (ssum != tsum) return {false, COST(0)};
        if (Dinic(sg, s, t) < ssum) return {false, COST(0)};

        // come down to min-cost circulation
        G.init(V);
        for (int i = 0; i < (int)edges.size(); i++) {
            auto &e = edges[i];
            const auto &ge = sg.get_edge(i);
            G.add_edge(ge.from, ge.to, ge.cap, ge.flow, e.cost);
        }
        MinCostCirculation(G);

        // find min-cost
        COST res = 0;
        for (int i = 0; i < (int)edges.size(); i++) {
            auto &e = edges[i];
            const auto &ge = G.get_edge(i);
            e.flow = e.upper_cap - ge.cap;
            res += e.flow * e.cost;
        }

        // find dual
        dual.resize(V);
        while (true) {
            bool update = false;
            for (int v = 0; v < V; v++) {
                for (const auto &e : G[v]) {
                    if (!e.cap) continue;
                    auto cost2 = dual[v] + e.cost;
                    if (cost2 < dual[e.to]) {
                        dual[e.to] = cost2;
                        update = true;
                    }
                }
            }
            if (!update) break;
        }
        return {true, res};
    }
};


//------------------------------//
// Solver
//------------------------------//

// Yosupo Libray Checker - Minimum Cost b-flow
void Yosupo_Minimum_Cost_b_flow() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);

    int N, M;
    cin >> N >> M;
    MinCostBFlow<long long, i128> mcf(N);
    vector<i128> B(N);
    for (int i = 0; i < N; i++) cin >> B[i], mcf.add_ds(i, B[i]);
    vector<int> s(M), t(M);
    vector<i128> l(M), u(M), c(M);
    for (int i = 0; i < M; i++) {
        cin >> s[i] >> t[i] >> l[i] >> u[i] >> c[i];
        mcf.add_edge(s[i], t[i], l[i], u[i], c[i]);
    }

    auto [exist, res] = mcf.solve();
    if (!exist) cout << "infeasible" << '\n';
    else {
        cout << res << '\n';
        for (int i = 0; i < N; i++) cout << mcf.dual[i] << '\n';
        for (int i = 0; i < M; i++) cout << mcf.edges[i].flow << '\n';
    }
}


int main() {
    Yosupo_Minimum_Cost_b_flow();
}