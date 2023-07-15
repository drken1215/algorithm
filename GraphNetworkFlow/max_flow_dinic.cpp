//
// max-flow (Dinic's algorithm)
//
// verified
//   AOJ Course GRL_6_A Network Flow - Maximum Flow
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_6_A&lang=jp
//
//   AtCoder Library Practice Contest D - Maxflow
//     https://atcoder.jp/contests/practice2/tasks/practice2_d
//
//   ABC 259 G - Grid Card Game
//     https://atcoder.jp/contests/abc259/tasks/abc259_g
//


#include <bits/stdc++.h>
using namespace std;


// edge class (for network-flow)
template<class FLOWTYPE> struct FlowEdge {
    // core members
    int rev, from, to;
    FLOWTYPE cap, icap, flow;
    
    // constructor
    FlowEdge(int r, int f, int t, FLOWTYPE c)
    : rev(r), from(f), to(t), cap(c), icap(c), flow(0) {}
    void reset() { cap = icap, flow = 0; }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowEdge& E) {
        return s << E.from << "->" << E.to << '(' << E.flow << '/' << E.icap << ')';
    }
};

// graph class (for network-flow)
template<class FLOWTYPE> struct FlowGraph {
    // core members
    vector<vector<FlowEdge<FLOWTYPE>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    
    // constructor
    FlowGraph(int n = 0) : list(n) { }
    void init(int n = 0) {
        list.assign(n, FlowEdge<FLOWTYPE>());
        pos.clear();
    }
    
    // getter
    vector<FlowEdge<FLOWTYPE>> &operator [] (int i) {
        return list[i];
    }
    const vector<FlowEdge<FLOWTYPE>> &operator [] (int i) const {
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowEdge<FLOWTYPE> &get_rev_edge(const FlowEdge<FLOWTYPE> &e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    FlowEdge<FLOWTYPE> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowEdge<FLOWTYPE> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowEdge<FLOWTYPE>> get_edges() const {
        vector<FlowEdge<FLOWTYPE>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowEdge<FLOWTYPE> &e : list[i]) e.reset();
        }
    }
    void change_edge(FlowEdge<FLOWTYPE> &e, FLOWTYPE new_cap, FLOWTYPE new_flow) {
        FlowEdge<FLOWTYPE> &re = get_rev_edge(e);
        e.cap = new_cap - new_flow, e.icap = new_cap, e.flow = new_flow;
        re.cap = new_flow;
    }
    
    // add_edge
    void add_edge(int from, int to, FLOWTYPE cap) {
        pos.emplace_back(from, (int)list[from].size());
        list[from].push_back(FlowEdge<FLOWTYPE>((int)list[to].size(), from, to, cap));
        list[to].push_back(FlowEdge<FLOWTYPE>((int)list[from].size() - 1, to, from, 0));
    }

    // debug
    friend ostream& operator << (ostream& s, const FlowGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

template<class FLOWTYPE> FLOWTYPE Dinic
 (FlowGraph<FLOWTYPE> &G, int s, int t, FLOWTYPE limit_flow)
{
    FLOWTYPE current_flow = 0;
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
            for (const FlowEdge<FLOWTYPE> &e : G[v]) {
                if (level[e.to] < 0 && e.cap > 0) {
                    level[e.to] = level[v] + 1;
                    if (e.to == t) return;
                    que.push(e.to);
                }
            }
        }
    };
    
    // Dinic DFS
    auto dfs = [&](auto self, int v, FLOWTYPE up_flow) {
        if (v == t) return up_flow;
        FLOWTYPE res_flow = 0;
        for (int &i = iter[v]; i < (int)G[v].size(); ++i) {
            FlowEdge<FLOWTYPE> &e = G[v][i], &re = G.get_rev_edge(e);
            if (level[v] >= level[e.to] || e.cap == 0) continue;
            FLOWTYPE flow = self(self, e.to, min(up_flow - res_flow, e.cap));
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
            FLOWTYPE flow = dfs(dfs, s, limit_flow - current_flow);
            if (!flow) break;
            current_flow += flow;
        }
    }
    return current_flow;
};

template<class FLOWTYPE> FLOWTYPE Dinic(FlowGraph<FLOWTYPE> &G, int s, int t) {
    return Dinic(G, s, t, numeric_limits<FLOWTYPE>::max());
}



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// AOJ
void AOJ_Course_GRL_6_A() {
    int V, E;
    cin >> V >> E;
    FlowGraph<int> G(V);
    for (int i = 0; i < E; ++i) {
        int u, v, c;
        cin >> u >> v >> c;
        G.add_edge(u, v, c);
    }
    int res = Dinic(G, 0, V-1);
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


int main() {
    //AOJ_Course_GRL_6_A();
    ACL_practice_D();
    //ABC_259_G();
}

