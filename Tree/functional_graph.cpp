//
// (連結とは限らない) Functional Graph をサイクルと森に分解する
//
// verified:
//   AtCoder ABC 256 E - Takahashi's Anguish
//     https://atcoder.jp/contests/abc256/tasks/abc256_e
//
//   AtCoder ABC 387 F - Count Arrays
//     https://atcoder.jp/contests/abc387/tasks/abc387_f
//


#include <bits/stdc++.h>
using namespace std;


// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
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

// 連結とは限らない Functional Graph を、サイクルと森に分解していく
// G[v] := 頂点 v から出ている辺
template<class T = long long> struct RunFunctionalGraph {
    // cycles
    const int NOT_IN_CYCLE = -1;
    vector<vector<int>> roots;  // nodes in each cycle
    vector<vector<Edge<T>>> cycles;  // the cycles
    vector<int> cmp;  // order in tye cycle
    
    // trees
    vector<vector<Edge<T>>> childs;
    vector<unordered_map<int,int>> id;  // id[v][w] := the index of node w in G[v]
    vector<long long> siz;  // the size of v-subtree
    
    // for finding lca
    vector<vector<int>> parent;
    vector<int> root, depth, which_cycle;

    // Euler tour
    vector<int> tour; // the node-number of i-th element of Euler-tour
    vector<int> v_s_id, v_t_id; // the index of Euler-tour of node v
    vector<int> e_id; // the index of edge e (v*2 + (0: root to leaf, 1: leaf to root)

    // constructor
    RunFunctionalGraph() { }
    RunFunctionalGraph(const Graph<T> &G) {
        init(G);
    }

    // get first / last id of node v in Euler tour
    int vs(int v) { return v_s_id[v]; }
    int vt(int v) { return v_t_id[v]; }
    int get_v(int id) { return tour[id]; }

    // get edge-id of (pv, v) in Euler tour
    int e(int v, bool leaf_to_root = false) {
        assert(cmp[v] == NOT_IN_CYCLE);
        if (!leaf_to_root) return e_id[v * 2];
        else return e_id[v * 2 + 1];
    }
    int e(int u, int v) {
        if (depth[u] < depth[v]) return e(v);
        else return e(u, false);
    }
    pair<int, int> get_e(int id) { 
        return make_pair(tour[id], tour[id + 1]);
    }

    // get_parent(v, p) := the parent of v directed for p
    int get_parent(int v) { return parent[0][v];  }
    int get_root(int v) { return root[v]; }
    int kth_ancestor(int v, int k) {
        if (k > depth[v]) return root[v];
        int goal_depth = depth[v] - k;
        for (int i = (int)parent.size()-1; i >= 0; i--)
            if (parent[i][v] != -1 && depth[parent[i][v]] >= goal_depth) 
                v = parent[i][v];
        return v;
    }
    int get_parent(int v, int p) {
        assert(v != p && root[v] == root[p]);
        int lca = get_lca(v, p);
        if (lca != v) return parent[0][v];
        else return kth_ancestor(p, depth[p] - depth[v] - 1);
    }

    // lca(u, v)
    int get_lca(int u, int v) {
        assert(root[u] == root[v]);
        if (depth[u] > depth[v]) swap(u, v);
        for (int i = 0; i < (int)parent.size(); i++) {
            if ((depth[v] - depth[u]) & (1<<i))
                v = parent[i][v];
        }
        if (u == v) return u;
        for (int i = (int)parent.size()-1; i >= 0; i--) {
            if (parent[i][u] != parent[i][v]) {
                u = parent[i][u];
                v = parent[i][v];
            }
        }
        return parent[0][u];
    }

    // dist(u, v)
    long long get_dist(int u, int v) {
        assert(which_cycle[u] == which_cycle[v]);
        if (root[u] == root[v]) {
            int lca = get_lca(u, v);
            return depth[u] + depth[v] - depth[lca] * 2;
        } else {
            int res = depth[u] + depth[v];
            u = root[u], v = root[v];
            int cycledis = max(cmp[u], cmp[v]) - min(cmp[u], cmp[v]);
            cycledis = min(cycledis, (int)cycles[which_cycle[v]].size() - cycledis);
            res += cycledis;
            return res;
        }
    }

    // is node v in s-t path?
    bool is_on_path(int s, int t, int v) {
        return get_dist(s, v) + get_dist(v, t) == get_dist(s, t);
    };
    
    // init
    void detect_all_cycles(const Graph<T> &G) {
        int N = (int)G.size();
        roots.clear(), cycles.clear();
        cmp.assign(N, NOT_IN_CYCLE);
        vector<bool> seen(N, false), finished(N, false);
        vector<int> history;
        auto pop_history = [&]() -> void {
            while (!history.empty()) {
                int v = history.back();
                finished[v] = true;
                history.pop_back();
            }
        };
        auto detect_a_node_in_the_cycle = [&](int v) {
            do {
                seen[v] = true;
                history.push_back(v);
                v = G[v][0].to;
                if (finished[v]) {
                    v = -1;
                    break;
                }
            } while (!seen[v]);
            pop_history();
            return v;
        };
        auto reconstruct = [&](int r) -> pair<vector<int>, vector<Edge<T>>> {
            vector<int> sub_roots;
            vector<Edge<T>> cycle;
            int v = r, iter = 0;
            do {
                sub_roots.emplace_back(v);
                cycle.emplace_back(G[v][0]);
                cmp[v] = iter++;
                v = G[v][0].to;
            } while (v != r);
            return {sub_roots, cycle};
        };
        for (int v = 0; v < (int)G.size(); ++v) {
            if (finished[v]) continue;
            int r = detect_a_node_in_the_cycle(v);
            if (r == -1) continue;
            auto [sub_roots, cycle] = reconstruct(r);
            if (!cycle.empty()) {
                roots.emplace_back(sub_roots);
                cycles.emplace_back(cycle);
            }
        }
    }
    void init(const Graph<T> &G, int s = 0) {
        int N = (int)G.size();

        // step 0: assertion
        for (int v = 0; v < N; v++) assert(G[v].size() == 1);
        
        // step 1: detect all cycles
        detect_all_cycles(G);

        // step 2: construct trees
        // min non-negative i such that n <= 2^i
        int D = 0;
        while ((1LL << D) < N) D++;
        parent.assign(D + 1, vector<int>(N, -1)), childs.resize(N);
        for (int v = 0; v < N; v++) {
            if (cmp[v] != NOT_IN_CYCLE) {
                parent[0][v] = v;
            } else {
                childs[G[v][0].to].emplace_back(Edge<T>(G[v][0].to, v, G[v][0].val));
                parent[0][v] = G[v][0].to;
            }
        }
        for (int i = 0; i < D; i++) for (int v = 0; v < N; v++) {
            parent[i + 1][v] = parent[i][parent[i][v]];
        }

        // step 3: run trees
        depth.resize(N), root.resize(N), which_cycle.resize(N);
        siz.resize(N), id.resize(N);
        tour.resize(N * 2 - 1), v_s_id.resize(N), v_t_id.resize(N), e_id.resize(N * 2);
        int ord = 0;
        auto rec = [&](auto &&rec, int v, int d, int cid, int r) -> int {
            int sum = 1;
            depth[v] = d, root[v] = r, which_cycle[v] = cid;
            ord++;
            for (int i = 0; i < (int)childs[v].size(); i++) {
                int ch = childs[v][i].to;
                id[v][ch] = i;
                e_id[ch * 2] = ord - 1;
                sum += rec(rec, ch, d + 1, cid, r);
                tour[ord] = v, v_t_id[v] = ord, e_id[ch * 2 + 1] = ord - 1;
                ord++;
            }
            siz[v] = sum;
            return sum;
        };
        for (int cid = 0; cid < (int)roots.size(); cid++) {
            for (auto r : roots[cid]) {
                ord = 0;
                rec(rec, r, 0, cid, r);
            }
        }
    }
};


//------------------------------//
// Examples
//------------------------------//

// AtCoder ABC 256 E - Takahashi's Anguish
void ABC_256_E() {
    int N; cin >> N;
    long long res = 0;
    Graph<long long> G(N);
    vector<long long> X(N), C(N);
    for (int i = 0; i < N; i++) cin >> X[i], X[i]--;
    for (int i = 0; i < N; i++) cin >> C[i];
    for (int i = 0; i < N; i++) G.add_edge(i, X[i], C[i]);
    RunFunctionalGraph<long long> fg(G);
    for (const auto &cycle : fg.cycles) {
        long long sub = 1LL << 60;
        for (const auto &e : cycle) sub = min(sub, e.val);
        res += sub;
    }
    cout << res << endl;
}

// AtCoder ABC 387 F - Count Arrays
template<class T_VAL, class T_MOD>
constexpr T_VAL mod_inv(T_VAL a, T_MOD m) {
    T_VAL b = m, u = 1, v = 0;
    while (b > 0) {
        T_VAL t = a / b;
        a -= t * b, swap(a, b);
        u -= t * v, swap(u, v);
    }
    u %= m;
    if (u < 0) u += m;
    return u;
}
template<int MOD = 998244353, bool PRIME = true> struct Fp {
    // inner value
    unsigned int val;
    
    // constructor
    constexpr Fp() : val(0) { }
    template<std::signed_integral T> constexpr Fp(T v) {
        long long tmp = (long long)(v % (long long)(get_umod()));
        if (tmp < 0) tmp += get_umod();
        val = (unsigned int)(tmp);
    }
    template<std::unsigned_integral T> constexpr Fp(T v) {
        val = (unsigned int)(v % get_umod());
    }
    constexpr long long get() const { return val; }
    constexpr static int get_mod() { return MOD; }
    constexpr static unsigned int get_umod() { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator + () const { return Fp(*this); }
    constexpr Fp operator - () const { return Fp() - Fp(*this); }
    constexpr Fp operator + (const Fp &r) const { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp &r) {
        val += r.val;
        if (val >= get_umod()) val -= get_umod();
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) {
        val -= r.val;
        if (val >= get_umod()) val += get_umod();
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) {
        unsigned long long tmp = val;
        tmp *= r.val;
        val = (unsigned int)(tmp % get_umod());
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) {
        return *this = *this * r.inv(); 
    }
    constexpr Fp pow(long long n) const {
        assert(n >= 0);
        Fp res(1), mul(*this);
        while (n) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const {
        if (PRIME) {
            assert(val);
            return pow(get_umod() - 2);
        } else {
            assert(val);
            return mod_inv((long long)(val), get_umod());
        }
    }

    // other operators
    constexpr bool operator == (const Fp &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const {
        return this->val != r.val;
    }
    constexpr bool operator < (const Fp &r) const {
        return this->val < r.val;
    }
    constexpr bool operator > (const Fp &r) const {
        return this->val > r.val;
    }
    constexpr bool operator <= (const Fp &r) const {
        return this->val <= r.val;
    }
    constexpr bool operator >= (const Fp &r) const {
        return this->val >= r.val;
    }
    constexpr Fp& operator ++ () {
        ++val;
        if (val == get_umod()) val = 0;
        return *this;
    }
    constexpr Fp& operator -- () {
        if (val == 0) val = get_umod();
        --val;
        return *this;
    }
    constexpr Fp operator ++ (int) {
        Fp res = *this;
        ++*this;
        return res;
    }
    constexpr Fp operator -- (int) {
        Fp res = *this;
        --*this;
        return res;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD> &x) {
        long long tmp = 1;
        is >> tmp;
        tmp = tmp % (long long)(get_umod());
        if (tmp < 0) tmp += get_umod();
        x.val = (unsigned int)(tmp);
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD> &x) {
        return os << x.val;
    }
    friend constexpr Fp<MOD> pow(const Fp<MOD> &r, long long n) {
        return r.pow(n);
    }
    friend constexpr Fp<MOD> inv(const Fp<MOD> &r) {
        return r.inv();
    }
};
void ABC_387_F() {
    using mint = Fp<>;
    int N, M;
    cin >> N >> M;
    vector<long long> A(N);
    Graph<long long> G(N);
    for (int i = 0; i < N; i++) {
        cin >> A[i], A[i]--;
        G.add_edge(i, A[i]);
    }
    RunFunctionalGraph<long long> fg(G);

    vector dp(N, vector(M, mint(0))), sdp(N, vector(M+1, mint(0)));
    auto rec = [&](auto &&rec, int v) -> void {    
        for (auto e : fg.childs[v]) rec(rec, e.to);
        for (long long val = 0; val < M; val++) {
            dp[v][val] = 1;
            for (auto e : fg.childs[v]) dp[v][val] *= sdp[e.to][val+1];
            sdp[v][val+1] = sdp[v][val] + dp[v][val];
        }
    };
    mint res = 1;
    for (auto rts : fg.roots) {
        mint sub = 0;
        for (auto r : rts) rec(rec, r);
        for (long long val = 0; val < M; val++) {
            mint tmp = 1;
            for (auto r : rts) tmp *= dp[r][val];
            sub += tmp;
        }
        res *= sub;
    }
    cout << res << endl;
}


int main() {
    //ABC_256_E();
    ABC_387_F();
}