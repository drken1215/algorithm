//
// min-cut (Dinic's algorithm)
//
// example:
//   KUPC 2016 E - 柵
//     https://atcoder.jp/contests/kupc2016/tasks/kupc2016_e
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

// Dinic
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

void KUPC_2016_E() {
    const vector<int> dx = {1, 0, -1, 0};
    const vector<int> dy = {0, 1, 0, -1};
    
    int H, W;
    cin >> H >> W;
    vector<string> S(H);
    for (int i = 0; i < H; ++i) cin >> S[i];
    
    // (i,j) を表す index
    auto ind = [&](int i, int j) -> int { return i * W + j; };
    
    // フローネットワークの構築
    const int INF = 1<<29;
    int N = H * W;
    int s = N * 2, t = N * 2 + 1;
    FlowGraph<int> G(N * 2 + 2);
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            // in -> out
            G.add_edge(ind(i, j), ind(i, j) + N, (S[i][j] == '.' ? 1 : INF));
            
            // s -> ヤギ
            if (S[i][j] == 'X') {
                G.add_edge(s, ind(i, j), INF);
            }
            
            // 外周 -> t
            if (i == 0 || i == H-1 || j == 0 || j == W-1) {
                G.add_edge(ind(i, j) + N, t, INF);
            }
            
            // 上下左右
            for (int d = 0; d < 4; ++d) {
                int i2 = i + dx[d], j2 = j + dy[d];
                if (i2 < 0 || i2 >= H || j2 < 0 || j2 >= W) continue;
                G.add_edge(ind(i, j) + N, ind(i2, j2), INF);
            }
        }
    }
    
    int res = Dinic(G, s, t);
    cout << (res < N ? res : -1) << endl;
}


int main() {
    KUPC_2016_E();
}


