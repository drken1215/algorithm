//
// 最小費用 b-flow by ネットワーク単体法
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

// Network Simplex Method
template<class FLOW, class COST> struct NetworkSimplex {
    struct Edge {
        int from, to;
        FLOW cap;
        COST cost;
        friend ostream& operator << (ostream& s, const Edge& e) {
            return s << e.from << " -> " << e.to << " (" << e.cap << ", " << e.cost << ")";
        }
    };
    struct Parent {
        int p, e;
        FLOW up, down;
    };

    // inner values
    int N, M;
    vector<Edge> edges;
    vector<FLOW> lowers;  // for edge i
    vector<FLOW> dss;  // for node v (demand and supply)

    // intermediate results
    int BUCKET_SIZE, MINOR_LIMIT;
    vector<Parent> parents;
    vector<int> depth, nex, pre, candidates;

    // results
    bool feasible;
    COST total_cost;
    vector<FLOW> flows;
    vector<COST> pots;

    // constructor
    explicit NetworkSimplex(int n = 0) : N(n), dss(n) {}

    // debugger
    friend ostream& operator << (ostream&s, const NetworkSimplex &ns) {
        s << "feasibility: " << (ns.feasible ? "true" : "false") << '\n';
        s << "optimal value: " << ns.total_cost << '\n';
        for (int i = 0; i < ns.M; i += 2) {
            auto e = ns.edges[i];
            s << e.from << " -> " << e.to << ": " << ns.flows[i / 2]
              << " (remained cap: " << e.cap << ", cost: " << e.cost << ")\n";
        }
        for (int v = 0; v < ns.N; v++) cout << "node " << v << ": " << ns.pots[v] << '\n';
        return s;
    }

    // setter
    void add_edge(int from, int to, FLOW lower, FLOW upper, COST cost) {
        edges.push_back({from, to, upper - lower, cost});
        edges.push_back({to, from, 0, -cost});
        lowers.push_back(lower);
        dss[from] -= lower;
        dss[to] += lower;
        M = (int)edges.size();
    }
    void add_ds(int v, FLOW ds) {
        assert(v >= 0 && v < N);
        dss[v] += ds;
    }

    // solver
    pair<bool, COST> solve() {
        BUCKET_SIZE = max(int(sqrt(double(M)) * 0.2), 10);
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
                if (edges[ei].cap) {
                    COST clen = edges[ei].cost + pots[edges[ei ^ 1].to] - pots[edges[ei].to];
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
        pots.assign(N + 1, 0); 
        parents.resize(N), depth.assign(N + 1, 1); 
        nex.assign((N + 1) * 2, 0), pre.assign((N + 1) * 2, 0);
        COST inf_cost = 1;
        for (int i = 0; i < M; i += 2) {
            inf_cost += (edges[i].cost >= 0 ? edges[i].cost : -edges[i].cost);
        }
        edges.reserve(M + N * 2);
        for (int i = 0; i < N; i++) {
            if (dss[i] >= 0) {
                edges.push_back({i, N, 0, inf_cost});
                edges.push_back({N, i, dss[i], -inf_cost});
                pots[i] = -inf_cost;
            } else {
                edges.push_back({i, N, -dss[i], -inf_cost});
                edges.push_back({N, i, 0, inf_cost});
                pots[i] = inf_cost;
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
            if (dss[i] >= 0) {
                if (edges[M + i * 2 + 1].cap) feasible = false;
            } else {
                if (edges[M + i * 2].cap) feasible = false;
            }
        }
        if (!feasible) return false;
        total_cost = 0;
        flows.clear();
        for (int i = 0; i < M; i += 2) {
            flows.push_back(lowers[i / 2] + edges[i ^ 1].cap);
            total_cost += flows.back() * edges[i].cost;
        }
        pots.pop_back();
        return true;
    }

    void push_flow(int ei0) {
        int u0 = edges[ei0 ^ 1].to, v0 = edges[ei0].to, del_u = v0;
        FLOW flow = edges[ei0].cap;
        COST clen = edges[ei0].cost + pots[u0] - pots[v0];
        bool del_u_side = true;
        int lca = get_lca(u0, v0, flow, del_u_side, del_u);
        if (flow) {
            int u = u0, v = v0;
            while (u != lca) parents[u].up += flow, parents[u].down -= flow, u = parents[u].p;
            while (v != lca) parents[v].up -= flow, parents[v].down += flow, v = parents[v].p;
        }
        int u = u0, par = v0;
        auto p_caps = make_pair(edges[ei0].cap - flow, edges[ei0 ^ 1].cap + flow);
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
                if (idx % 2 == 0) d++, pots[idx / 2] += p_diff, depth[idx / 2] = d;
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
            if (!edges[ei].cap) {
                swap(candidates[i], candidates.back());
                candidates.pop_back();
                continue;
            }
            COST clen = edges[ei].cost + pots[edges[ei ^ 1].to] - pots[edges[ei].to];
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


//------------------------------//
// Solver
//------------------------------//

// Yosupo Libray Checker - Minimum Cost b-flow
void Yosupo_Minimum_Cost_b_flow() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);

    int N, M;
    cin >> N >> M;
    NetworkSimplex<long long, i128> ns(N);
    vector<i128> B(N);
    for (int i = 0; i < N; i++) cin >> B[i], ns.add_ds(i, B[i]);
    vector<int> s(M), t(M);
    vector<i128> l(M), u(M), c(M);
    for (int i = 0; i < M; i++) {
        cin >> s[i] >> t[i] >> l[i] >> u[i] >> c[i];
        ns.add_edge(s[i], t[i], l[i], u[i], c[i]);
    }

    auto [exist, res] = ns.solve();
    if (!exist) cout << "infeasible" << '\n';
    else {
        cout << res << '\n';
        for (int i = 0; i < N; i++) cout << ns.pots[i] << '\n';
        for (int i = 0; i < M; i++) cout << ns.flows[i] << '\n';
    }
}


int main() {
    Yosupo_Minimum_Cost_b_flow();
}