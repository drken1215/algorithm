//
// 最小費用 b-flow by network simplex method
//
// verified
//   Yosupo Library Checker - Minimum Cost b-flow (N <= 100)
//   (primal-dual: 負閉路 NG, cost-scaling: 15 ms, simplex: 2 ms)
//     https://judge.yosupo.jp/problem/min_cost_b_flow
//
//   ABC 421 G - Increase to make it Increasing (N <= 300)
//   (primal-dual: 4 ms, cost-scaling: 24 ms, simplex: 2 ms)
//     https://atcoder.jp/contests/abc421/tasks/abc421_g
//
//   KUPC 2014 I - Rain (N <= 10^4, ただし K <= 15 によって流量が小さい)
//   (primal-dual: 35 ms, cost-scaling: TLE, simplex: 681 ms)
//     https://atcoder.jp/contests/kupc2014/tasks/kupc2014_i
//
//   JAG 夏合宿 2013 Day4 I - Multi Path Story (N <= 1000)
//   (primal-dual: 18 ms, cost-scaling: 2590 ms, simplex: 6 ms)
//     https://onlinejudge.u-aizu.ac.jp/problems/2627
//
//   Educational Codeforces Round 80 F. Red-Blue Graph (N <= 200)
//   (primal-dual: 62 ms, cost-scaling: 46 ms, simplex: 62 ms)
//     https://codeforces.com/contest/1288/problem/F
//
//   UTPC 2011 H - キャッシュ戦略 (N <= 10^4, 本来は流量 M <= 10 だがそれを活かさない解法をしている）
//   (primal-dual: TLE, cost-scaling: TLE, simplex: 1130 ms)
//   #pragma GCC optimize("Ofast") を入れるとむしろ遅くなる！！！
//     https://atcoder.jp/contests/utpc2011/tasks/utpc2011_8
//
//   Codeforces Round 826 (Div. 3) G. Kirill and Company 
//   (primal-dual: 負閉路 NG, cost-scaling: 1046 ms, simplex: 62 ms)
//     https://codeforces.com/contest/1741/problem/G
//
//   AtCoder ABC 393 G - Unevenness (for using frac<i128> and dual)
//   (primal-dual: 負閉路 NG, cost-scaling: 830 ms, simplex: 116 ms)
//     https://atcoder.jp/contests/abc393/tasks/abc393_g
//


#include <bits/stdc++.h>
using namespace std;


// output stream
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class S, class T> ostream& operator << (ostream &s, const pair<S, T> &P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 2> &P)
{ return s << '<' << P[0] << "," << P[1] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 3> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 4> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << "," << P[3] << '>'; }
template<class T> ostream& operator << (ostream &s, const vector<string> &P)
{ for (int i = 0; i < P.size(); ++i) { s << P[i] << endl; } return s; }
template<class T> ostream& operator << (ostream &s, const vector<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const deque<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const vector<vector<T>> &P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, const set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const multiset<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const unordered_set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const unordered_map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }


//--------------------------------//
// b-flow
//--------------------------------//

// Network Simplex Method
template<class FLOW, class COST> struct MinCostBFlowByNetworkSimplex {
    // inner Edge
    struct InnerEdge {
        int from, to;
        FLOW cap;
        COST cost;
        InnerEdge(int from_, int to_, FLOW cap_, COST cost_) : from(from_), to(to_), cap(cap_), cost(cost_) {}
    };
    struct Parent {
        int p, e;
        FLOW up, down;
    };

    // inner values
    int N, original_edge_size;
    vector<InnerEdge> edges;
    vector<FLOW> dss;  // demand (< 0) and supply (> 0)
    bool feasible;
    COST total_cost;
    vector<COST> dual;

    // intermediate results
    int BUCKET_SIZE, MINOR_LIMIT;
    vector<Parent> parents;
    vector<int> depth, nex, pre, candidates;

    // constructor
    explicit MinCostBFlowByNetworkSimplex(int n = 0) : N(n), dss(n) {}

    // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        edges.emplace_back(from, to, cap, cost);
        edges.emplace_back(to, from, 0, -cost);
    }
    void set_ds(int v, FLOW ds) {
        assert(0 <= v && v < N);
        dss[v] = ds;
    }
    void set_ds(const vector<FLOW> &vds) {
        assert((int)vds.size() == N);
        dss = vds;
    }

    // getter
    FLOW get_flow(int i) const {
        return edges[(i * 2) ^ 1].cap;
    }
    COST get_dual(int v) const {
        return dual[v];
    }
    vector<COST> get_duals() const {
        return dual;
    }

    // solver
    pair<bool, COST> solve() {
        BUCKET_SIZE = max(int(sqrt(double(edges.size())) * 0.2), 10);
        MINOR_LIMIT = max(int(BUCKET_SIZE * 0.1), 3);
        precompute();
        candidates.reserve(BUCKET_SIZE);
        int ei = 0;
        while (true) {
            for (int i = 0; i < MINOR_LIMIT; i++) if (!minor()) break;
            COST best = 0;
            int best_ei = -1;
            candidates.clear();
            for (int i = 0; i < (int)edges.size(); i++) {
                if (edges[ei].cap > 0) {
                    COST clen = edges[ei].cost + dual[edges[ei ^ 1].to] - dual[edges[ei].to];
                    if (clen < 0) {
                        if (clen < best) best = clen, best_ei = ei;
                        candidates.push_back(ei);
                        if ((int)candidates.size() == BUCKET_SIZE) break;
                    }
                }
                ei++;
                if (ei == (int)edges.size()) ei = 0;
            }
            if (candidates.empty()) break;
            push_flow(best_ei);
        }
        if (!postcompute()) return {false, COST(-1)};
        else return {true, total_cost};
    }

    void connect(int a, int b) {
        nex[a] = b, pre[b] = a;
    }

    void precompute() {
        original_edge_size = (int)edges.size();
        dual.assign(N + 1, 0); 
        parents.resize(N), depth.assign(N + 1, 1); 
        nex.assign((N + 1) * 2, 0), pre.assign((N + 1) * 2, 0);
        COST inf_cost = 1;
        for (int i = 0; i < (int)edges.size(); i += 2) {
            inf_cost += (edges[i].cost >= 0 ? edges[i].cost : -edges[i].cost);
        }
        edges.reserve((int)edges.size() + N * 2);
        for (int i = 0; i < N; i++) {
            if (dss[i] >= 0) {
                edges.push_back(InnerEdge(i, N, 0, inf_cost));
                edges.push_back(InnerEdge(N, i, dss[i], -inf_cost));
                dual[i] = -inf_cost;
            } else {
                edges.push_back(InnerEdge(i, N, -dss[i], -inf_cost));
                edges.push_back(InnerEdge(N, i, 0, inf_cost));
                dual[i] = inf_cost;
            }
            int e = (int)edges.size() - 2;
            parents[i] = {N, e, edges[e].cap, edges[e ^ 1].cap};
        }
        depth[N] = 0;
        for (int i = 0; i < N + 1; i++) connect(i * 2, i * 2 + 1);
        for (int i = 0; i < N; i++) connect(i * 2 + 1, nex[N * 2]), connect(N * 2, i * 2);
    }

    bool postcompute() {
        for (int i = 0; i < N; i++) {
            edges[parents[i].e].cap = parents[i].up;
            edges[parents[i].e ^ 1].cap = parents[i].down;
        }
        feasible = true;
        for (int i = 0; i < N; i++) {
            int e = original_edge_size + i * 2;
            if (dss[i] >= 0) {
                if (edges[e ^ 1].cap > 0) feasible = false;
            } else {
                if (edges[e].cap > 0) feasible = false;
            }
        }
        if (!feasible) return false;
        total_cost = 0;
        for (int i = 0; i < (int)edges.size(); i += 2) {
            total_cost += edges[i ^ 1].cap * edges[i].cost;
        }
        dual.pop_back();
        return true;
    }

    void push_flow(int ei0) {
        int u0 = edges[ei0 ^ 1].to, v0 = edges[ei0].to, del_u = v0;
        FLOW f = edges[ei0].cap;
        COST clen = edges[ei0].cost + dual[u0] - dual[v0];
        bool del_u_side = true;
        int lca = get_lca(u0, v0, f, del_u_side, del_u);
        if (f > 0) {
            int u = u0, v = v0;
            while (u != lca) parents[u].up += f, parents[u].down -= f, u = parents[u].p;
            while (v != lca) parents[v].up -= f, parents[v].down += f, v = parents[v].p;
        }
        int u = u0, par = v0;
        auto p_caps = make_pair(edges[ei0].cap - f, edges[ei0 ^ 1].cap + f);
        COST p_diff = -clen;
        if (!del_u_side) {
            swap(u, par); 
            swap(p_caps.first, p_caps.second);
            p_diff *= -1;
        }
        int par_e = ei0 ^ (del_u_side ? 0 : 1);
        while (par != del_u) {
            int d = depth[par], idx = u * 2;
            while (idx != u * 2 + 1) {
                if (idx % 2 == 0) d++, dual[idx / 2] += p_diff, depth[idx / 2] = d;
                else d--;
                idx = nex[idx];
            }
            connect(pre[u * 2], nex[u * 2 + 1]);
            connect(u * 2 + 1, nex[par * 2]);
            connect(par * 2, u * 2);
            swap(parents[u].e, par_e);
            par_e ^= 1;
            swap(parents[u].up, p_caps.first); 
            swap(parents[u].down, p_caps.second);
            swap(p_caps.first, p_caps.second);
            int next_u = parents[u].p; 
            parents[u].p = par;
            par = u;
            u = next_u;
        }
        edges[par_e].cap = p_caps.first;
        edges[par_e ^ 1].cap = p_caps.second;
    }

    bool minor() {
        if (candidates.empty()) return false;
        COST best = 0;
        int best_ei = -1;
        int i = 0;
        while (i < int(candidates.size())) {
            int ei = candidates[i];
            if (edges[ei].cap <= 0) {
                swap(candidates[i], candidates.back());
                candidates.pop_back();
                continue;
            }
            COST clen = edges[ei].cost + dual[edges[ei ^ 1].to] - dual[edges[ei].to];
            if (clen >= 0) {
                swap(candidates[i], candidates.back());
                candidates.pop_back();
                continue;
            }
            if (clen < best) best = clen, best_ei = ei;
            i++;
        }
        if (best_ei == -1) return false;
        push_flow(best_ei);
        return true;
    }

    int get_lca(int u, int v, FLOW &flow, bool &del_u_side, int &del_u) {
        auto up_u = [&]() {
            if (parents[u].down < flow) flow = parents[u].down, del_u = u, del_u_side = true;
            u = parents[u].p;
        };
        auto up_v = [&]() {
            if (parents[v].up <= flow) flow = parents[v].up, del_u = v, del_u_side = false;
            v = parents[v].p;
        };
        if (depth[u] >= depth[v]) {
            int num = depth[u] - depth[v];
            for (int i = 0; i < num; i++) up_u();
        } else {
            int num = depth[v] - depth[u];
            for (int i = 0; i < num; i++) up_v();
        }
        while (u != v) up_u(), up_v();
        return u;
    }
};

// b-flow manager
template<class FLOW, class COST> struct MinCostBFlow {
    // Edge
    struct InnerEdge {
        int from, to;
        FLOW lower_cap, upper_cap, flow;
        COST cost;
        InnerEdge(int from_, int to_, FLOW lower_, FLOW upper_, COST cost_)
            : from(from_), to(to_), lower_cap(lower_), upper_cap(upper_), flow(0), cost(cost_) {}
        friend ostream& operator << (ostream& s, const InnerEdge& e) {
            return s << e.from << "->" << e.to 
            << " (" << e.flow << "/" << e.lower_cap << "~" << e.upper_cap << ", " << e.cost << ")";
        }
    };

    // inner values
    int N;
    vector<InnerEdge> edges;
    vector<FLOW> lower_dss, upper_dss, dss;  // demand (< 0) and supply (> 0)
    vector<COST> dual;
    
    // constructor
    explicit MinCostBFlow(int n = 0) : N(n), lower_dss(n, 0), upper_dss(n, 0), dss(n, 0) {}

    // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        edges.push_back(InnerEdge(from, to, 0, cap, cost));
    }
    void add_edge(int from, int to, FLOW lower_cap, FLOW upper_cap, COST cost) {
        assert(lower_cap <= upper_cap);
        edges.push_back(InnerEdge(from, to, lower_cap, upper_cap, cost));
    }
    void set_ds(int v, FLOW ds) {
        assert(0 <= v && v < N);
        lower_dss[v] = ds, upper_dss[v] = ds;
    }
    void set_ds(int v, FLOW lower_ds, FLOW upper_ds) {
        assert(0 <= v && v < N);
        assert(lower_ds <= upper_ds);
        lower_dss[v] = lower_ds, upper_dss[v] = upper_ds;
    }

    // getter
    InnerEdge &get_edge(int i) {
        return edges[i];
    }
    const InnerEdge &get_edge(int i) const {
        return edges[i];
    }
    vector<InnerEdge> get_edges() const {
        return edges;
    }
    COST get_dual(int v) const {
        return dual[v];
    }
    vector<COST> get_duals() const {
        return dual;
    }

    // solver
    bool pre_compute() {
        bool need_super_node = false;
        for (int v = 0; v < N; v++) {
            if (lower_dss[v] == upper_dss[v]) dss[v] = lower_dss[v];
            else need_super_node = true;
        }

        // lower_ds, upper_ds -> strict ds
        if (need_super_node) {
            int super = N;
            dss.assign(N + 1, 0);
            for (int v = 0; v < N; v++) {
                if (lower_dss[v] >= 0) {
                    add_edge(super, v, lower_dss[v], upper_dss[v], 0);
                } else if (upper_dss[v] < 0) {
                    add_edge(v, super, -upper_dss[v], -lower_dss[v], 0);
                } else {
                    add_edge(super, v, upper_dss[v], 0);
                    add_edge(v, super, -lower_dss[v], 0);
                }
            }
        }

        // push lower_cap
        for (const auto &e : edges) {
            dss[e.to] += e.lower_cap, dss[e.from] -= e.lower_cap;
        }
        return need_super_node;
    }
    pair<bool, COST> solve(const string solver = "network_simplex", bool calc_potential = false) {
        bool need_super_node = pre_compute();
        COST res = 0;
        if (solver == "network_simplex") {
            MinCostBFlowByNetworkSimplex<FLOW, COST> G(N + (int)need_super_node);
            G.set_ds(dss);
            for (const auto &e : edges) G.add_edge(e.from, e.to, e.upper_cap - e.lower_cap, e.cost);
            auto [feasible, mincost] = G.solve();
            if (!feasible) return {false, COST(0)};
            for (int i = 0; i < (int)edges.size(); i++) {
                auto &e = edges[i];
                e.flow = e.lower_cap + G.get_flow(i);
                res += e.flow * e.cost;
            }
            if (calc_potential) {
                dual = G.get_duals();
                if (need_super_node) dual.pop_back();
            }
        }
        return {true, res};
    }
};


//------------------------------//
// Solver
//------------------------------//

// Yosupo Libray Checker - Minimum Cost b-flow
using i128 = __int128_t;
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
void Yosupo_Minimum_Cost_b_flow(const string solver) {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);

    int N, M;
    cin >> N >> M;
    MinCostBFlow<long long, i128> G(N);
    vector<i128> B(N);
    for (int i = 0; i < N; i++) cin >> B[i], G.set_ds(i, B[i]);
    vector<int> s(M), t(M);
    vector<i128> l(M), u(M), c(M);
    for (int i = 0; i < M; i++) {
        cin >> s[i] >> t[i] >> l[i] >> u[i] >> c[i];
        G.add_edge(s[i], t[i], l[i], u[i], c[i]);
    }
    auto [exist, res] = G.solve(solver, true);
    if (!exist) cout << "infeasible" << '\n';
    else {
        const auto &dual = G.get_duals();
        const auto &es = G.get_edges();
        cout << res << '\n';
        for (auto v : dual) cout << v << '\n';
        for (const auto &e : es) cout << e.flow << '\n';
    }
}

// ABC 421 G - Increase to make it Increasing
void ABC_421_G(const string solver) {
    long long N, M, INF = 1LL<<45; cin >> N >> M;
    vector<long long> A(N), D(N+1, INF);
    D[0] = 0;
    for (int i = 0; i < N; i++) {
        cin >> A[i];
        if (i) D[i] = A[i] - A[i-1];
    }

    MinCostBFlow<long long, long long> G(N+1);
    for (int v = 0; v <= N; v++) {
        if (D[v] >= 0) G.set_ds(v, 0, D[v]);
        else G.set_ds(v, -INF, D[v]);
    }
    for (int i = 0; i < M; i++) {
        long long u, v; cin >> u >> v; u--;
        G.add_edge(v, u, INF, 1);
    }
    auto [flag, cost] = G.solve(solver);
    cout << (flag ? cost : -1) << endl;
}

// KUPC 2014 I - Rain
void KUPC_2014_I(const string solver) {
    int N, M, K, INF = 1LL<<20;
    cin >> N >> M >> K;
    vector<long long> C(K), A(M), B(M), D(M), b(N, 0); 
    for (int i = 0; i < K; i++) cin >> C[i], C[i]--;
    for (int i = 0; i < M; i++) cin >> A[i] >> B[i] >> D[i], A[i]--, B[i]--;
    for (int i = 0; i < K; i++) b[A[C[i]]]++, b[B[C[i]]]--;
    MinCostBFlow<int, long long> G(N);
    for (int i = 0; i < N; i++) G.set_ds(i, b[i]);
    for (int i = 0; i < M; i++) G.add_edge(B[i], A[i], INF, D[i]);
    auto [flag, cost] = G.solve(solver);
    cout << (flag ? cost : -1) << endl;
}

// JAG 夏合宿 2013 Day4 I - Multi Path Story (AOJ 2627)
void AOJ_2627(const string solver) {
    long long N, INF = 1LL << 45;
    cin >> N;
    MinCostBFlow<long long, long long> G(N + 1);
    long long t = N;
    for (int v = 0; v < N; v++) {
        int D;
        cin >> D;
        for (int i = 0; i < D; i++) {
            long long to, w;
            cin >> to >> w, to--;
            G.add_edge(v, to, 1, INF, w);
        }
    }
    for (int v = 1; v < N; v++) G.add_edge(v, t, 0, INF, 0);
    G.add_edge(t, 0, 0, INF, 0);
    auto [flag, cost] = G.solve(solver);
    cout << cost << endl;
}

// Educational Codeforces Round 80 F. Red-Blue Graph
void EducationalCodeforces80_F(const string solver) {
    long long L, R, M, costR, costB, INF = 10000;
    string sl, sr;
    cin >> L >> R >> M >> costR >> costB >> sl >> sr;
    MinCostBFlow<long long, long long> G(L + R);
    for (int i = 0; i < L; i++) {
        if (sl[i] == 'R') G.set_ds(i, 1, INF);
        else if (sl[i] == 'B') G.set_ds(i, -INF, -1);
        else G.set_ds(i, -INF, INF);
    }
    for (int j = 0; j < R; j++) {
        if (sr[j] == 'R') G.set_ds(j+L, -INF, -1);
        else if (sr[j] == 'B') G.set_ds(j+L, 1, INF);
        else G.set_ds(j+L, -INF, INF);
    }
    for (int i = 0; i < M; i++) {
        int u, v;
        cin >> u >> v, u--, v--;
        G.add_edge(u, v+L, 1, costR);
        G.add_edge(v+L, u, 1, costB);
    }
    auto [flag, mincost] = G.solve(solver);
     if (!flag) {
        cout << -1 << endl;
        return;
     }
    auto es = G.get_edges();
    string res(M, 'U');
    for (int i = 0; i < M; i++) {
        auto e = es[i*2], re = es[i*2+1];
        if (e.flow == 1) res[i] = 'R';
        else if (re.flow == 1) res[i] = 'B';
    }
    cout << mincost << endl;
    cout << res << endl;
}

// UTPC 2011 H - キャッシュ戦略
void UTPC2011_H(const string solver) {
    int M, N, K;
    cin >> M >> N >> K;
    vector<int> W(N), A(K);
    for (int i = 0; i < N; i++) cin >> W[i];
    for (int i = 0; i < K; i++) cin >> A[i], A[i]--;
    int s = K * 3, t = s + 1;
    MinCostBFlow<int, int> G(K * 3 + 2);
    for (int v = 0; v < K; v++) {
        G.add_edge(s, v, 1, W[A[v]]);
        G.add_edge(v+K*2, v, 1, W[A[v]]);  // 空の状態からは常に行ける
        G.add_edge(v+K, t, 1, 0);
        G.add_edge(v, v+K, 1, 1, 0);  // 流量下限も 1

        // 次の同じ色のボールが来るまでキープする場合
        for (int v2 = v+1; v2 < K; v2++) {
            if (A[v] == A[v2]) {
                G.add_edge(v+K, v2, 1, 0);
                break;
            }
        }

        // 箱からボールを取り出す頂点との絡み
        if (v+1 < K) {
            G.add_edge(v+K, (v+1)+K*2, 1, 0);
            G.add_edge(v+K*2, (v+1)+K*2, M, 0);
        }
    }
    G.add_edge(t, s, M, 0);
    auto [flag, mincost] = G.solve(solver);
    cout << mincost << endl;
}

// Codeforces Round 826 (Div. 3) G. Kirill and Company
void Codeforces826_G(const string solver) {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int T;
    cin >> T;
    while (T--) {
        int N, M, F, K, INF = 6;
        cin >> N >> M;
        vector<vector<int>> G(N);
        for (int i = 0; i < M; i++) {
            int a, b; cin >> a >> b; a--, b--;
            G[a].emplace_back(b), G[b].emplace_back(a);
        }
        vector<vector<int>> prev(N);
        vector<int> dp(N, -1);
        queue<int> que;
        dp[0] = 0;
        que.push(0);
        while (!que.empty()) {
            auto v = que.front(); que.pop();
            for (auto v2 : G[v]) {
                if (dp[v2] == -1) {
                    dp[v2] = dp[v] + 1;
                    prev[v2].emplace_back(v);
                    que.push(v2);
                } else if (dp[v2] == dp[v] + 1) {
                    prev[v2].emplace_back(v);
                }
            }
        }
        cin >> F;
        vector<int> where(F), allnum(N, 0), nothavenum(N, 0), havenum(N, 0);
        for (int i = 0; i < F; i++) {
            cin >> where[i], where[i]--;
            allnum[where[i]]++;
        }
        cin >> K;
        for (int i = 0; i < K; i++) {
            int v; cin >> v, v--;
            nothavenum[where[v]]++;
        }
        for (int v = 0; v < N; v++) havenum[v] = allnum[v] - nothavenum[v];

        int s = N * 2, t = N;
        MinCostBFlow<int, int> FG(N * 2 + 1);
        for (int v = 0; v < N; v++) {
            if (havenum[v] > 0) FG.add_edge(s, v, havenum[v], 0);
            if (nothavenum[v] > 0) FG.add_edge(v, v+N, 1, -nothavenum[v]);
            FG.add_edge(v, v+N, INF, 0);
            for (auto v2 : prev[v]) FG.add_edge(v+N, v2, INF, 0);
        }
        FG.add_edge(t, s, INF, 0);
        auto [flag, mincost] = FG.solve(solver);
        int res = K + mincost;
        cout << res << endl;
    }
}

// AtCoder ABC 393 G - Unevenness
template<class T> struct SternBrocotTree {
    template<class Func> static tuple<T, T, T, T> binary_search(Func check, T lim) {
        assert(check(0, 1));
        assert(!check(1, 0));
        auto rec = [&](auto &&rec, bool which, T &a, T &b, T c, T d) -> void {
            if (a + c > lim || b + d > lim) return;
            if (check(a + c, b + d) == which) {
                a += c, b += d;
                rec(rec, which, a, b, c + c, d + d);
            }
            if (a + c <= lim && b + d <= lim && check(a + c, b + d) == which) a += c, b += d;
        };
        T a = 0, b = 1, c = 1, d = 0;
        while (a + c <= lim && b + d <= lim) {
            rec(rec, true, a, b, c, d);
            rec(rec, false, c, d, a, b);
        }
        return {a, b, c, d};
    }
};
template<class T = long long> struct frac {
    // gcd
    static T gcd(T a, T b) {
        a = max(a, -a), b = max(b, -b);
        while (b) {
            a %= b;
            swap(a, b);
        }
        return a;
    }

    // inner values
    T first, second;

    // constructor
    frac& normalize() {
        if (first == 0 && second != 0) {
            second = 1;
            return *this;
        }
        if (second == 0 && first != 0) {
            first = 1;
            return *this;
        }
        if (second < 0) first = -first, second = -second;
        T d = gcd(max(first, -first), second);
        if (d == 0) first = 0, second = 1;
        else first /= d, second /= d;
        return *this;
    }
    frac(const frac&) = default;
    frac& operator = (const frac&) = default;
    constexpr frac(T f = 0, T s = 1) : first(f), second(s) { 
        normalize(); 
    }
    constexpr frac& operator = (T a) { 
        *this = frac(a, 1); 
        return *this;
    }
    constexpr long double to_double() const {
        assert(second != 0);
        return (long double)(first) / (long double)(second);
    }
    friend constexpr long double to_double(const frac &r) {
        return r.to_double();
    }

    // comparison operators
    constexpr bool operator == (const frac &r) const {
        return this->first == r.first && this->second == r.second;
    }
    constexpr bool operator != (const frac &r) const {
        return this->first != r.first || this->second != r.second;
    }
    constexpr bool operator < (const frac &r) const {
        return this->first * r.second < this->second * r.first;
    }
    constexpr bool operator > (const frac &r) const {
        return this->first * r.second > this->second * r.first;
    }
    constexpr bool operator <= (const frac &r) const {
        return this->first * r.second <= this->second * r.first;
    }
    constexpr bool operator >= (const frac &r) const {
        return this->first * r.second >= this->second * r.first;
    }
    
    // arithmetic operators
    constexpr frac& operator += (const frac &r) {
        this->first = this->first * r.second + this->second * r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator -= (const frac &r) {
        this->first = this->first * r.second - this->second * r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator *= (const frac &r) {
        this->first *= r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator /= (const frac &r) {
        this->first *= r.second;
        this->second *= r.first;
        this->normalize();
        return *this;
    }
    constexpr frac operator + () const { return frac(*this); }
    constexpr frac operator - () const { return frac(0) - frac(*this); }
    constexpr frac operator + (const frac &r) const { return frac(*this) += r; }
    constexpr frac operator - (const frac &r) const { return frac(*this) -= r; }
    constexpr frac operator * (const frac &r) const { return frac(*this) *= r; }
    constexpr frac operator / (const frac &r) const { return frac(*this) /= r; }
    friend constexpr ostream& operator << (ostream &os, const frac<T> &x) {
        os << x.first; 
        if (x.second != 1) os << "/" << x.second;
        return os;
    }
};
void ABC_393_G(const string solver) {
    using i128 = __int128_t;
    using FR = frac<i128>;
    using SBT = SternBrocotTree<long long>;
    long long N, P, Q;
    cin >> N >> P >> Q;
    FR K(P, Q);
    vector A(N, vector(N, 0LL));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) cin >> A[i][j];

    // 目的関数の値
    auto calc_obj = [&](const vector<FR> &x) -> FR {
        FR res = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            if (i+1 < N) res += max(x[i*N+j] - x[(i+1)*N+j], x[(i+1)*N+j] - x[i*N+j]);
            if (j+1 < N) res += max(x[i*N+j] - x[i*N+j+1], x[i*N+j+1] - x[i*N+j]);
        }
        return res;
    };

    // 小さい λ では負の値になり、ある程度大きい λ では 0 になる。0 になる瞬間が最適解。
    auto calc_penalty = [&](const vector<FR> &x) -> FR {
        FR res = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            res += max(x[i*N+j] - FR(A[i][j]), FR(A[i][j]) - x[i*N+j]);
        }
        return res - K;
    };

    // Stern-Brocot 木上の二分探索を実施する
    /*
    min: Σ_{u < v}(max(0, x[v] - x[u]) + max(0, x[u] - x[v]))
            + λ(Σ_{v}(max(0, x[v] - A[v]) + max(0, A[v] - x[v])) - P/Q)
    */
    auto optimize = [&](FR r, vector<FR> &x) -> FR {
        MinCostBFlow<FR, FR> G(N * N + 1);
        int s = N * N;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            if (i+1 < N) {
                int u = i*N+j, v = (i+1)*N+j;
                G.add_edge(u, v, FR(1), FR(0));
                G.add_edge(v, u, FR(1), FR(0));
            }
            if (j+1 < N) {
                int u = i*N+j, v = i*N+j+1;
                G.add_edge(u, v, FR(1), FR(0));
                G.add_edge(v, u, FR(1), FR(0));
            }
            int v = i*N+j;
            G.add_edge(s, v, r, FR(A[i][j]));
            G.add_edge(v, s, r, -FR(A[i][j]));
        }
        auto [flag, cost] = G.solve(solver, true);  // dual も求める
        auto ans = G.dual;
        assert(ans.size() == N * N + 1);
        for (int i = 0; i < N * N; i++) x[i] = ans[i] - ans[s];
        return cost - r * K;
    };
    auto check = [&](i128 a, i128 b) -> bool {
        FR r(a, b);
        vector<FR> x(N*N);
        if (a == 1 && b == 0) return false;
        auto cost = optimize(r, x);
        return calc_penalty(x) > 0;
    };
    
    vector<FR> pre_x(N*N), nex_x(N*N), x(N*N);
    if (!check(0, 1)) {
        vector<FR> x(N*N);
        optimize(FR(0, 1), x);
    } else {
        auto [a, b, c, d] = SBT::binary_search(check, 10000000000000LL);
        FR pre_r(a, b), nex_r(c, d);
        auto pre_all_cost = optimize(pre_r, pre_x);
        auto nex_all_cost = optimize(nex_r, nex_x);
        auto pre_obj = calc_obj(pre_x), nex_obj = calc_obj(nex_x);
        auto pre_penalty = calc_penalty(pre_x), nex_penalty = calc_penalty(nex_x);
        for (int i = 0; i < N*N; i++) {
            x[i] = (pre_x[i]*(-nex_penalty) + nex_x[i]*pre_penalty) / (pre_penalty - nex_penalty);
        }
    }
    auto res = calc_obj(x);
    cout << fixed << setprecision(20) << to_double(res) << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) cout << to_double(x[i*N+j]) << " ";
        cout << endl;
    }
}


int main() {
    //Yosupo_Minimum_Cost_b_flow("network_simplex");
    //ABC_421_G("network_simplex");
    //KUPC_2014_I("network_simplex");
    //AOJ_2627("network_simplex");
    //EducationalCodeforces80_F("network_simplex");
    //UTPC2011_H("network_simplex");
    //Codeforces826_G("network_simplex");
    ABC_393_G("network_simplex");
}