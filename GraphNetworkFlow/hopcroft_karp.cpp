//
// Hopcroft-Karp の最大二部マッチング
//
// verified
//   AOJ 1163 カードゲーム
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1163&lang=jp
//
//   TTPC 2022 H - Colorful Graph
//     https://atcoder.jp/contests/ttpc2022/tasks/ttpc2022_h
//


#include <bits/stdc++.h>
using namespace std;


struct HopcroftKarp {
    // input
    int size_left, size_right;
    vector<vector<int>> list; // left to right

    // results
    vector<int> lr, rl;
    
    // intermediate results
    vector<bool> seen, matched;
    vector<int> level;
    
    // constructor
    HopcroftKarp(int l, int r) : size_left(l), size_right(r), list(l, vector<int>()) { }
    void add_edge(int from, int to) {
        list[from].push_back(to);
    }

    // getter, debugger
    const vector<int> &operator [] (int i) const { 
        return list[i];
    }
    friend ostream& operator << (ostream& s, const HopcroftKarp& G) {
        s << endl;
        for (int i = 0; i < G.list.size(); ++i) {
            s << i << ": ";
            for (int j = 0; j < G.list[i].size(); ++j) {
                s << G.list[i][j];
                if (j + 1 != G.list[i].size()) s << ", ";
            }
            s << endl;
        }
        return s;
    }
    
    // solver
    void hobfs() {
        queue<int> que;
        for (int left = 0; left < size_left; ++left) {
            level[left] = -1;
            if (!matched[left]) {
                que.push(left);
                level[left] = 0;
            }
        }
        level[size_left] = size_left;
        while (!que.empty()) {
            int left = que.front();
            que.pop();
            for (int i = 0; i < list[left].size(); ++i) {
                int right = list[left][i];
                int next = rl[right];
                if (level[next] == -1) {
                    level[next] = level[left] + 1;
                    que.push(next);
                }
            }
        }
    }
    bool hodfs(int left) {
        if (left == size_left) return true;
        if (seen[left]) return false;
        seen[left] = true;
        for (int i = 0; i < list[left].size(); ++i) {
            int right = list[left][i];
            int next = rl[right];
            if (next == -1) next = size_left;
            if (level[next] > level[left] && hodfs(next)) {
                rl[right] = left;
                return true;
            }
        }
        return false;
    }
    int solve() {
        seen.assign(size_left, false);
        matched.assign(size_left, false);
        level.assign(size_left + 1, -1);
        lr.assign(size_left, -1);
        rl.assign(size_right, -1);
        int res = 0;
        while (true) {
            hobfs();
            seen.assign(size_left, false);
            bool finished = true;
            for (int left = 0; left < size_left; ++left) {
                if (!matched[left] && hodfs(left)) {
                    matched[left] = true;
                    ++res;
                    finished = false;
                }
            }
            if (finished) break;
        }
        for (int r = 0; r < size_right; r++) {
            if (rl[r] != -1) lr[rl[r]] = r;
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

// AOJ 1163 カードゲーム
void AOJ_1163() {
    int N, M;
    while (cin >> N >> M, N) {
        HopcroftKarp G(N, M);
        vector<int> left(N), right(M);
        for (int i = 0; i < N; ++i) cin >> left[i];
        for (int i = 0; i < M; ++i) cin >> right[i];
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                if (gcd(left[i], right[j]) > 1) G.add_edge(i, j);
            }
        }
        cout << G.solve() << endl;
    }
}

// TTPC 2022 H - Colorful Graph
// Edge Class
template<class T> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1), val(-1) { }
    Edge(int f, int t, T v = -1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& E) {
        return s << E.from << "->" << E.to;
    }
};

// graph class
template<class T> struct Graph {
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    
    Graph(int n = 0) : list(n), reversed_list(n) { }
    void init(int n = 0) {
        list.assign(n, vector<Edge<T>>());
        reversed_list.assign(n, vector<Edge<T>>());
    }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    const vector<Edge<T>> &get_rev_edges(int i) const { return reversed_list[i]; }
    const size_t size() const { return list.size(); }
        
    void add_edge(int from, int to, T val = -1) {
        list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }
    
    void add_bidirected_edge(int from, int to, T val = -1) {
        list[from].push_back(Edge(from, to, val));
        list[to].push_back(Edge(to, from, val));
        reversed_list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }

    friend ostream &operator << (ostream &s, const Graph &G) {
        s << endl;
        for (int i = 0; i < G.size(); ++i) {
            s << i << " -> ";
            for (const auto &e : G[i]) s << e.to << " ";
            s << endl;
        }
        return s;
    }
};

// strongly connected components decomposition
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
void TTPC_2022_H() {
    int N, M, a, b;
    cin >> N >> M;
    Graph<int> G(N);
    for (int i = 0; i < M; i++) {
        cin >> a >> b;
        a--, b--;
        G.add_edge(a, b, i);
    }
    SCC scc(G);

    const auto &dag = scc.dag;
    const auto &groups = scc.groups;

    int V = dag.size();
    HopcroftKarp hk(V, V);
    for (int v = 0; v < V; v++) {
        vector<bool> can(V, false);
        can[v] = true;
        queue<int> que;
        que.push(v);
        while (!que.empty()) {
            auto x = que.front();
            que.pop();
            for (auto e : dag[x]) {
                if (can[e.to]) continue;
                can[e.to] = true;
                que.push(e.to);
            }
        }
        for (int v2 = 0; v2 < V; v2++) {
            if (v2 != v && can[v2]) hk.add_edge(v, v2);
        }
    }
    int max_flow = hk.solve();

    vector<int> source, ans(V, -1);
    for (int v = 0; v < V; v++) if (hk.rl[v] == -1) source.push_back(v);
    for (int c = 0; c < source.size(); c++) {
        int v = source[c];
        while (v != -1) {
            ans[v] = c + 1;
            v = hk.lr[v];
        }
    }

    vector<int> res(N, -1);
    for (int v = 0; v < V; v++) for (auto v2 : groups[v]) res[v2] = ans[v];
    for (int i = 0; i < N; i++) cout << res[i] << " ";
    cout << endl;
}


int main() {
    AOJ_1163();
    //TTPC_2022_H();
}