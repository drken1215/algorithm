//
// max-flow (Dinic's algorithm)
//
// verified
//   典型アルゴリズム問題集 上級〜エキスパート編 E - 最大流
//     https://atcoder.jp/contests/pastbook2022/tasks/pastbook2022_e
//     https://atcoder.jp/contests/tessoku-book/tasks/tessoku_book_bp
//
//   AtCoder Library Practice Contest D - Maxflow
//     https://atcoder.jp/contests/practice2/tasks/practice2_d
//
//   ABC 259 G - Grid Card Game
//     https://atcoder.jp/contests/abc259/tasks/abc259_g
//
//   JAG 夏合宿 2011 Day4 D - Box Witch (AOJ 2313) (for change_edge)
//     https://onlinejudge.u-aizu.ac.jp/problems/2313
//
//   code festival 2014 上海 D - Maze (for reconstruct)
//     https://atcoder.jp/contests/code-festival-2014-china-open/tasks/code_festival_china_d
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


// edge class (for max-flow)
template<class FLOW> struct FlowEdge {
    // core members
    int rev, from, to;
    FLOW cap, icap, flow;
    
    // constructor
    constexpr FlowEdge() noexcept = default;
    constexpr FlowEdge(int rev, int from, int to, FLOW cap, FLOW rcap = 0) 
        : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(rcap) {
    }
    void reset() { 
        flow -= icap - cap;
        cap = icap;
    }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowEdge& e) {
        return s << e.from << " -> " << e.to << " (" << e.cap << ", " << e.flow << ")";
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
    void clear() {
        list.clear(), pos.clear();
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
    size_t size() const noexcept {
        return list.size();
    }
    FlowEdge<FLOW> &get_rev_edge(const FlowEdge<FLOW> &e) {
        return list[e.to][e.rev];
    }
    const FlowEdge<FLOW> &get_rev_edge(const FlowEdge<FLOW> &e) const {
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
    void reset() const {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowEdge<FLOW> &e : list[i]) e.reset();
        }
    }
    void change_edge(FlowEdge<FLOW> &e, FLOW new_cap, FLOW new_rcap) {
        assert(new_cap >= 0 && new_rcap >= 0);
        FlowEdge<FLOW> &re = get_rev_edge(e);
        e.cap = new_cap, e.icap = new_cap + new_rcap, e.flow = new_rcap;
        re.cap = new_rcap, re.icap = new_cap + new_rcap, re.flow = new_cap;
    }
    
    // add_edge
    void add_edge(int from, int to, FLOW cap, FLOW rcap = 0) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowEdge<FLOW>(to_id, from, to, cap, rcap));
        list[to].push_back(FlowEdge<FLOW>(from_id, to, from, rcap, cap));
    }
    void add_bidirected_edge(int from, int to, FLOW cap) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        add_edge(from, to, cap, cap);
    }

    // augment
    FLOW augment(int s, int t, FLOW up_flow = numeric_limits<FLOW>::max()) {
        vector<bool> seen(size(), false);
        auto dfs = [&](auto &&dfs, int v, FLOW up_flow) -> FLOW {
            if (v == t) return up_flow;
            seen[v] = true;
            for (int i = 0; i < (int)list[v].size(); i++) {
                FlowEdge<FLOW> &e = list[v][i], &re = get_rev_edge(e);
                if (seen[e.to] || e.cap <= 0) continue;
                FLOW flow = dfs(dfs, e.to, min(up_flow, e.cap));
                if (flow > 0) {
                    e.cap -= flow, e.flow += flow;
                    re.cap += flow, re.flow -= flow;
                    return flow;
                }
            }  
            return FLOW(0); 
        };
        return dfs(dfs, s, up_flow);
    };

    // find reachable nodes from node s (1: s-domain, -1: t-domain, 0: no reach)
    vector<int> find_cut(int s, int t) const {
        vector<int> res(size(), 0);
        auto dfs_s = [&](auto &&dfs_s, int v) -> void {
            res[v] = 1;
            for (const auto &e : list[v]) {
                if (res[e.to] || e.cap <= 0) continue;
                dfs_s(dfs_s, e.to);
            }
        };
        auto dfs_t = [&](auto &&dfs_t, int v) -> void {
            res[v] = -1;
            for (const auto &e : list[v]) {
                auto re = get_rev_edge(e);
                if (res[e.to] || re.cap <= 0) continue;
                dfs_t(dfs_t, e.to);
            }
        };
        dfs_s(dfs_s, s), dfs_t(dfs_t, t);
        return res;
    }

    // check if the s-t flow is feasible
    bool is_feasible(int s, int t) const {
        vector<FLOW> b(list.size(), FLOW(0));
        for (int v = 0; v < (int)list.size(); v++) {
            for (const auto &e : list[v]) {
                b[v] += (e.flow - get_rev_edge(e).flow) / 2;
            }
        }
        if (b[s] + b[t] != 0) return false;
        for (int v = 0; v < (int)list.size(); v++) {
            if (v != s && v != t && b[v] != FLOW(0)) return false;
        }
        return true;
    }
    bool is_feasible(int s, int t, FLOW flow) const {
        vector<FLOW> b(list.size(), FLOW(0));
        for (int v = 0; v < (int)list.size(); v++) {
            for (const auto &e : list[v]) {
                b[v] += (e.flow - get_rev_edge(e).flow) / 2;
            }
        }
        if (b[s] != flow) return false;
        if (b[t] != -flow) return false;
        for (int v = 0; v < (int)list.size(); v++) {
            if (v != s && v != t && b[v] != FLOW(0)) return false;
        }
        return true;
    }

    // decompose flow into s-t simple paths and cycles
    using Path = vector<FlowEdge<FLOW>>;
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
        vector<int> ptr(list.size(), 0), onpath(list.size(), -1);
        vector<pair<int, int>> route;
        vector<int> used;
        vector<Path> paths, cycles;

        auto next_arc = [&](int v) -> int {
            while (ptr[v] < (int)fg[v].size() && fg[v][ptr[v]].rem <= 0) ptr[v]++;
            return (ptr[v] < (int)fg[v].size() ? ptr[v] : -1);
        };
        auto extract = [&](int begin, bool is_cycle) {
            FLOW mi = numeric_limits<FLOW>::max();
            for (int k = begin; k < (int)route.size(); k++) {
                auto [v, i] = route[k];
                mi = min(mi, fg[v][i].rem);
            }
            vector<FlowEdge<FLOW>> seq;
            for (int k = begin; k < (int)route.size(); k++) {
                auto [v, i] = route[k];
                fg[v][i].rem -= mi;
                FlowEdge<FLOW> e = list[v][fg[v][i].eidx];
                e.flow = mi;
                seq.push_back(e);
            }
            if (is_cycle) cycles.push_back(move(seq));
            else paths.push_back(move(seq));
        };
        auto walk = [&](int start, bool stop_at_t) {
            route.clear();
            int v = start;
            onpath[v] = 0;
            used.push_back(v);
            while (true) {
                int i = next_arc(v), u = fg[v][i].to;
                route.push_back({v, i});
                if (stop_at_t && u == t) {
                    extract(0, false);
                    break;
                }
                if (onpath[u] != -1) {
                    extract(onpath[u], true);
                    break;
                }
                onpath[u] = (int)route.size();
                used.push_back(u);
                v = u;
            }
            for (int w : used) onpath[w] = -1;
            used.clear();
        };

        // extract all s-t paths
        while (next_arc(s) != -1) walk(s, true);

        // decompose remained circulation into cycles
        for (int v = 0; v < (int)list.size(); v++) while (next_arc(v) != -1) walk(v, false);

        return {paths, cycles};
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
            if (level[v] >= level[e.to] || e.cap <= 0) continue;
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
            if (flow <= 0) break;
            current_flow += flow;
        }
    }
    return current_flow;
};

template<class FLOW> FLOW Dinic(FlowGraph<FLOW> &G, int s, int t) {
    return Dinic(G, s, t, numeric_limits<FLOW>::max());
}


//------------------------------//
// Examples
//------------------------------//

// 典型アルゴリズム問題集 上級〜エキスパート編 E - 最大流
void PAST_Max_Flow() {
    int V, E;
    cin >> V >> E;
    int s = 0, t = V - 1;
    FlowGraph<long long> G(V);
    for (int i = 0; i < E; ++i) {
        long long u, v, c;
        cin >> u >> v >> c, u--, v--;
        G.add_edge(u, v, c);
    }
    long long res = Dinic(G, s, t);

    /* debug: フローを復元した結果を示す */
    // auto [paths, cycles] = G.decompose(s, t);
    // for (int i = 0; i < (int)paths.size(); i++) {
    //     cout << "path " << i << ": " << paths[i][0].from;
    //     for (auto e : paths[i]) cout << " -> " << e.to;
    //     cout << " (" << paths[i][0].flow << ")" << endl;
    // }
    // for (int i = 0; i < (int)cycles.size(); i++) {
    //     cout << "cycle " << i << ": " << cycles[i][0].from;
    //     for (auto e : cycles[i]) cout << " -> " << e.to;
    //     cout << " (" << cycles[i][0].flow << ")" << endl;
    // }
    cout << res << endl;
}

// ACL practice D
void ACL_practice_D() {
    // 上下左右を表すベクトル
    const vector<int> DX = {1, 0, -1, 0};
    const vector<int> DY = {0, 1, 0, -1};
    
    // 入力受け取り
    int N, M;
    cin >> N >> M;
    vector<string> grid(N);
    for (int i = 0; i < N; ++i) cin >> grid[i];
    
    // フローネットワークを作る
    // 各マスの番号を 0, 1, ..., NM-1 とし、超頂点の番号を S = NM, T = NM+1 とする
    FlowGraph<int> G(N * M + 2);
    int S = N * M, T = N * M + 1;
    
    // マス (i, j) の頂点番号を返す関数
    auto index = [&](int i, int j) -> int { return i * M + j; };
    
    // 黒色マスと白色マスを結ぶ (黒色：i + j が偶数、白色：i + j が奇数)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            // 黒色マスならば、上下左右の 4 マスと辺を結んでいく
            if ((i + j) % 2 == 0 && grid[i][j] == '.') {
                for (int dir = 0; dir < 4; ++dir) {
                    int i2 = i + DX[dir], j2 = j + DY[dir];
                    if (i2 < 0 || i2 >= N || j2 < 0 || j2 >= M) continue;
                    
                    // どちらも空マスならば、ドミノを置けるので、辺を結ぶ
                    if (grid[i2][j2] == '.') {
                        G.add_edge(index(i, j), index(i2, j2), 1);
                    }
                }
            }
            
            // 超頂点 S から黒色マスへの辺を結ぶ
            if ((i + j) % 2 == 0 && grid[i][j] == '.') {
                G.add_edge(S, index(i, j), 1);
            }
            
            // 白色マスから超頂点 T への辺を結ぶ
            if ((i + j) % 2 == 1 && grid[i][j] == '.') {
                G.add_edge(index(i, j), T, 1);
            }
        }
    }
    
    // 最大流を流す
    int max_flow = Dinic(G, S, T);

    // フロー値が 1 となった辺を特定して、ドミノタイリングを復元する
    const auto &edges = G.get_edges();
    for (const auto &e : edges) {
        // 辺 e が超頂点に接続するものや、フロー値が 0 であるものはスキップ
        if (e.from == S || e.to == T || e.flow == 0) continue;
        
        // 辺 e の両端点に対応するマス
        int ifrom = e.from / M, jfrom = e.from % M;
        int ito = e.to / M, jto = e.to % M;
        
        // ドミノを置く
        if (ifrom == ito) {
            // ドミノを横に配置する場合
            if (jfrom > jto) swap(jfrom, jto);
            grid[ifrom][jfrom] = '>';
            grid[ito][jto] = '<';
        } else if (jfrom == jto) {
            // ドミノを縦に配置する場合
            if (ifrom > ito) swap(ifrom, ito);
            grid[ifrom][jfrom] = 'v';
            grid[ito][jto] = '^';
        }
    }
    
    // 出力
    cout << max_flow << endl;
    for (int i = 0; i < N; ++i) cout << grid[i] << endl;
}

// ABC 259 G
void ABC_259_G() {
    const long long INF = 1LL<<50;
    
    // 入力
    int H, W;
    cin >> H >> W;
    vector<vector<long long>> A(H, vector<long long>(W));
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) {
        cin >> A[i][j];
        A[i][j] = -A[i][j];
    }
    long long B = 0;
    vector<long long> S(H + W, 0);
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) S[i] += A[i][j];
    for (int j = 0; j < W; ++j) for (int i = 0; i < H; ++i) S[j+H] += A[i][j];
    for (int i = 0; i < H + W; ++i) B = min(B, S[i]);
    B = -B;
    
    // グラフを構築
    int source = H + W, sink = H + W + 1;
    FlowGraph<long long> G(H + W + 2);
    for (int i = 0; i < H; ++i) {
        G.add_edge(source, i, B);
        G.add_edge(i, sink, B + S[i]);
    }
    for (int j = 0; j < W; ++j) {
        G.add_edge(source, j+H, B + S[j+H]);
        G.add_edge(j+H, sink, B);
    }
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            long long cost = (A[i][j] <= 0 ? -A[i][j] : INF);
            G.add_edge(i, j+H, cost);
        }
    }
    long long flow = Dinic(G, source, sink);
    long long res = -(flow - B * (H + W));
    cout << res << endl;
}

// JAG 夏合宿 2011 Day4 D - Box Witch (AOJ 2313)
void AOJ_2313() {
    int N, M, Q, iter = 0;
    cin >> N >> M >> Q;
    vector<int> U(M), V(M), typ(Q), A(Q), B(Q);
    FlowGraph<int> G(N);
    map<pair<int,int>,int> ids;
    for (int i = 0; i < M; i++) {
        cin >> U[i] >> V[i], U[i]--, V[i]--;
        if (U[i] > V[i]) swap(U[i], V[i]);
        G.add_bidirected_edge(U[i], V[i], 1);
        ids[{U[i], V[i]}] = iter++;
    }
    for (int q = 0; q < Q; q++) {
        cin >> typ[q] >> A[q] >> B[q], A[q]--, B[q]--;
        if (A[q] > B[q]) swap(A[q], B[q]);
        if (!ids.count({A[q], B[q]})) {
            G.add_bidirected_edge(A[q], B[q], 0);
            ids[{A[q], B[q]}] = iter++;
        }
    }
    int s = 0, t = N-1;
    int flow = Dinic(G, s, t);
    for (int q = 0; q < Q; q++) {
        int eid = ids[{A[q], B[q]}];
        auto &e = G.get_edge(eid);
        assert(e.from == A[q] && e.to == B[q]); 
        if (typ[q] == 1) {
            assert(e.cap == 0 && G.get_rev_edge(e).cap == 0);
            G.change_edge(e, 1, 1);
        } else {
            assert(e.cap <= 2 && G.get_rev_edge(e).cap == 2 - e.cap);
            if (e.cap == 1) G.change_edge(e, 0, 0);
            else {
                // e が使われている場合を考える (閉路に含まれる場合と、s-t パスに含まれる場合がある)
                int from, to;
                if (e.cap == 0) from = e.from, to = e.to;
                else if (e.cap == 2) from = e.to, to = e.from;
                if (G.augment(from, to, 1) == 1) {
                    // e を含む閉路がある場合: その閉路を消せる
                    // ここで、e を含む s-t パスがある場合は、
                    // G.augment(from, to, 1) によって s-t パスが e を使わないものに張り変わる
                    G.change_edge(e, 0, 0);  // 最後に、e を消す
                } else {
                    // フローはパスと閉路に分解できることから、
                    // e を含む閉路がないならば、e が s-t パスに含まれることが保証される
                    // よって、残余グラフ上で t-s パスが存在することが保証される
                    // t から s へ逆向きに押し戻しておく (押し戻し時に e を戻すとは限らない)
                    assert(G.augment(t, s, 1) == 1);
                    flow--;
                    if (e.cap == 1) G.change_edge(e, 0, 0);
                    else {
                        // 今度は e を含む閉路の存在が保証されるので、上と同じことをする
                        if (e.cap == 0) from = e.from, to = e.to;
                        else if (e.cap == 2) from = e.to, to = e.from;
                        assert(G.augment(from, to, 1) == 1);
                        G.change_edge(e, 0, 0);
                    }
                }
            }
        }
        flow += G.augment(s, t);
        cout << flow << endl;
        assert(G.is_feasible(s, t, flow));
    }
}


int main() {
    PAST_Max_Flow();
    //ACL_practice_D();
    //ABC_259_G();
    //AOJ_2313(); 
}