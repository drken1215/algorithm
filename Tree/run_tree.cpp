//
// 木を走査して、さまざまな情報を求める
//
// verified:
//   Codeforces Round #614 (Div. 1) C. Xenon's Attack on the Gangs
//     https://codeforces.com/contest/1292/problem/C
//
//   CodeQUEEN 決勝 C - Path Intersection
//     https://atcoder.jp/contests/codequeen2023-final-open/tasks/codequeen2023_final_c
//
//   ABC 406 F - Compare Tree Weights (for Euler Tour)
//     https://atcoder.jp/contests/abc406/tasks/abc406_f
//


#include <bits/stdc++.h>
using namespace std;


// Run Tree (including Euler Tour)
template<class Graph = vector<vector<int>>> struct RunTree {
    // id[v][w] := the index of node w in G[v]
    vector<unordered_map<int, int>> id;

    // num[v][i] := the size of subtree of G[v][i] with parent v
    vector<vector<long long>> num;
    
    // for finding lca
    int root;
    vector<vector<int>> parent;
    vector<int> depth;

    // Euler tour
    vector<int> tour; // the node-number of i-th element of Euler-tour
    vector<int> v_s_id, v_t_id; // the index of Euler-tour of node v
    vector<int> e_id; // the index of edge e (v*2 + (0: root to leaf, 1: leaf to root))

    // constructor
    RunTree() {}
    RunTree(const Graph &G, int root = 0) : root(root) {
        init(G, root);
    }
    
    // init
    void init(const Graph &G, int root = 0) {
        int N = (int)G.size();
        id.assign(N, unordered_map<int,int>()), num.assign(N, vector<long long>());
        for (int v = 0; v < N; v++) num[v].assign((int)G[v].size(), 0);
        int h = 1, ord = 0;
        while ((1<<h) < N) h++;
        parent.assign(h, vector<int>(N, -1)), depth.resize(N);
        tour.resize(N*2-1), v_s_id.resize(N), v_t_id.resize(N), e_id.resize(N*2);
        rec(G, root, -1, 0, ord);
        for (int i = 0; i+1 < (int)parent.size(); ++i) {
            for (int v = 0; v < N; v++)
                if (parent[i][v] != -1)
                    parent[i+1][v] = parent[i][parent[i][v]];
        }
    }

    // get_size(u, v) := the size of subtree v with parent u
    long long get_size(int u, int v) {
        return num[u][id[u][v]];
    }

    // get first / last id of node v in Euler tour
    int vs(int v) { return v_s_id[v]; }
    int vt(int v) { return v_t_id[v]; }

    // get edge-id of (pv, v) in Euler tour
    int e(int v, bool leaf_to_root = false) {
        assert(v != root);
        if (!leaf_to_root) return e_id[v * 2];
        else return e_id[v * 2 + 1];
    }
    int e(int u, int v) {
        if (depth[u] < depth[v]) return e(v);
        else return e(u, false);
    }

    // lca(u, v)
    int get_lca(int u, int v) {
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
        int lca = get_lca(u, v);
        return depth[u] + depth[v] - depth[lca]*2;
    }

    // get_parent(v, p) := the parent of v directed for p
    int get_parent(int v, int p) {
        if (v == p) return -1;
        int lca = get_lca(v, p);
        if (lca != v) return parent[0][v];
        for (int i = (int)parent.size()-1; i >= 0; i--) {
            if (parent[i][p] != -1 && depth[parent[i][p]] > depth[v]) {
                p = parent[i][p];
            }
        }
        return p;
    }
    
    // rec
    int rec(const Graph &G, int v, int p, int d, int &ord) {
        int p_index = -1;
        int sum = 1;
        parent[0][v] = p, depth[v] = d;
        tour[ord] = v, v_s_id[v] = v_t_id[v] = ord;
        ord++;
        for (int i = 0; i < (int)G[v].size(); i++) {
            int ch = G[v][i];
            id[v][ch] = i;
            if (ch == p) {
                p_index = i;
                continue;
            }
            e_id[ch * 2] = ord - 1;
            int s = rec(G, ch, v, d+1, ord);
            num[v][i] = s;
            sum += s;
            tour[ord] = v;
            v_t_id[v] = ord;
            e_id[ch * 2 + 1] = ord - 1;
            ord++;
        }
        if (p_index != -1) num[v][p_index] = (int)G.size() - sum;
        return sum;
    }
};


//------------------------------//
// Examples
//------------------------------//

// Codeforces Round #614 (Div. 1) C. Xenon's Attack on the Gangs
void Codeforces_614_C() {
    int N;
    scanf("%d", &N);
    vector<vector<int>> G(N);
    for (int i = 0; i < N-1; ++i) {
        int u, v;
        scanf("%d %d", &u, &v);
        --u, --v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    
    // Run Tree
    RunTree rt(G);
    vector<vector<long long>> dp(N, vector<long long>(N, -1));
    
    // メモ化再帰
    auto rec = [&](auto self, int u, int v) -> long long {
        if (dp[u][v] != -1) return dp[u][v];
        if (dp[v][u] != -1) return dp[v][u];
        if (u == v) return 0;

        long long res = 0;
        int up = rt.get_parent(u, v), vp = rt.get_parent(v, u);
        res = max(res, self(self, up, v));
        res = max(res, self(self, vp, u));
        res += rt.get_size(up, u) * rt.get_size(vp, v);
        return dp[u][v] = dp[v][u] = res;
    };
    
    // 集計
    long long res = 0;
    for (int i = 0; i < N; ++i) for (int j = i+1; j < N; ++j) {
        res = max(res, rec(rec, i, j));
    }
    cout << res << endl;
}


// CodeQUEEN 決勝 C - Path Intersection
void CodeQUEEN_D() {
    int N, S, T;
    cin >> N >> S >> T;
    --S, --T;
    vector<vector<int>> G(N);
    for (int i = 0; i < N-1; ++i) {
        int u, v;
        cin >> u >> v;
        --u, --v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    
    RunTree rt(G, S);
    for (int v = 0; v < N; ++v) {
        int lca = rt.get_lca(v, T);
        cout << rt.get_dist(lca, v) + 1 << endl;
    }
}


// ABC 406 F - Compare Tree Weights
template <class Abel> struct BIT {
    Abel UNITY_SUM = 0;
    vector<Abel> dat;
    
    // [0, n)
    BIT(int n, Abel unity = 0) : UNITY_SUM(unity), dat(n, unity) { }
    void init(int n) {
        dat.assign(n, UNITY_SUM);
    }
    
    // a is 0-indexed
    inline void add(int a, Abel x) {
        for (int i = a; i < (int)dat.size(); i |= i + 1)
            dat[i] = dat[i] + x;
    }
    
    // [0, a), a is 0-indexed
    inline Abel sum(int a) {
        Abel res = UNITY_SUM;
        for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[i];
        return res;
    }
    
    // [a, b), a and b are 0-indexed
    inline Abel sum(int a, int b) {
        return sum(b) - sum(a);
    }
    
    // debug
    void print() {
        for (int i = 0; i < (int)dat.size(); ++i)
            cout << sum(i, i + 1) << ",";
        cout << endl;
    }
};
void ABC_406_F() {
    using Graph = vector<vector<int>>;
    int N, Q, u, v;
    cin >> N;
    Graph G(N);
    vector<pair<int,int>> edges(N-1);
    for (int i = 0; i < N-1; i++) {
        cin >> u >> v, u--, v--;
        edges[i] = {u, v};
        G[u].push_back(v), G[v].push_back(u);
    }
    RunTree rt(G);
    cin >> Q;
    long long all = 0;
    BIT<long long> bit(rt.tour.size() + 1);
    for (int qid = 0; qid < Q; qid++) {
        long long type, v, x, y;
        cin >> type;
        if (type == 1) {
            cin >> v >> x, v--;
            bit.add(rt.vs(v), x);
            all += x;
        } else if (type == 2) {
            cin >> y, y--;
            int u = edges[y].first, v = edges[y].second;
            if (rt.depth[u] > rt.depth[v]) swap(u, v);

            long long uv_size = rt.get_size(u, v), vu_size = rt.get_size(v, u);
            long long uv_sum = bit.sum(rt.vs(v), rt.vt(v)+1);  // +1 is necessary!
            long long vu_sum = all - uv_sum;
            long long uv = uv_size + uv_sum, vu = vu_size + vu_sum;
            long long res = abs(uv - vu);
            cout << res << '\n';
        } 
    }
}


int main() {
    //Codeforces_614_C();
    CodeQUEEN_D();
    //ABC_406_F();
}