//
// 2 変数劣モジュラ関数のグラフ表現
//
// verified:
//   AtCoder ARC 085 E - MUL
//     https://atcoder.jp/contests/arc085/tasks/arc085_c
//


/*
 N 個の 0-1 変数 x_0, x_1, ..., x_{N-1} について
 
 ・個別の変数 xi に関するコスト
    xi = 0 のとき false_cost, xi = 1 のとき true_cost
 
 ・2 変数 xi, xj 間の関係性についてのコスト
 　　(xi, xj) = (F, F): コスト A
 　　(xi, xj) = (F, T): コスト B
 　　(xi, xj) = (T, F): コスト C
 　　(xi, xj) = (T, T): コスト D
 　　(ただし、B + C >= A + D でなければならない)

 が設定されているときの最小コストを求める
 */


#include <bits/stdc++.h>
using namespace std;


// 2-variable submodular optimization
template<class T> struct TwoVariableSubmodularOpt {
    // edge class
    struct Edge {
        // core members
        int rev, from, to;
        T cap, icap, flow;
        
        // constructor
        Edge(int r, int f, int t, T c)
        : rev(r), from(f), to(t), cap(c), icap(c), flow(0) {}
        void reset() { cap = icap, flow = 0; }
        
        // debug
        friend ostream& operator << (ostream& s, const Edge& E) {
            return s << E.from << "->" << E.to << '(' << E.flow << '/' << E.icap << ')';
        }
    };
    
    // inner data
    int N, s, t;
    T offset;
    vector<vector<Edge>> list;
    vector<pair<int,int>> pos;
    
    // constructor
    TwoVariableSubmodularOpt() : N(0), s(0), t(0), offset(0) {}
    TwoVariableSubmodularOpt(int n) : N(n), s(n), t(n + 1), offset(0), list(n + 2) {}
    void init(int n = 0) {
        N = n, s = n, t = n + 1;
        offset = 0;
        list.assign(n + 2, Edge());
        pos.clear();
    }
    friend ostream& operator << (ostream& s, const TwoVariableSubmodularOpt &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }

    // add 1-Variable submodular functioin
    void add_single_cost(int xi, T false_cost, T true_cost) {
        assert(0 <= xi && xi < N);
        if (false_cost >= true_cost) {
            offset += true_cost;
            add_edge(s, xi, false_cost - true_cost);
        } else {
            offset += false_cost;
            add_edge(xi, t, true_cost - false_cost);
        }
    }
    
    // add constraint (xi = T, xj = F is penalty C)
    void add_penalty(int xi, int xj, T cost) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(cost >= 0);
        add_edge(xi, xj, cost);
    }
    
    // add constraint (xi = T, xj = T is only cost 0, and the others are cost C)
    void add_both_true_benefit(int xi, int xj, T cost) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(cost >= 0);
        add_edge(xj, 0, cost);
        add_edge(xi, xj, cost);
    }
    
    // add constraint (xi = F, xj = F is only cost 0, and the others are cost C)
    void add_both_false_benefit(int xi, int xj, T cost) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(cost >= 0);
        add_edge(xi, cost, 0);
        add_edge(xi, xj, cost);
    }
    
    // add general 2-value submodular function
    // (xi, xj) = (F, F): A, (F, T): B
    // (xi, xj) = (T, F): C, (T, T): D
    void add_submodular_function(int xi, int xj, T A, T B, T C, T D) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(B + C >= A + D);  // submodular constraint
        offset += A;
        add_single_cost(xi, 0, D - B);
        add_single_cost(xj, 0, B - A);
        add_penalty(xi, xj, B + C - A - D);
    }
    
    // solve
    T solve() {
        return dinic(s, t) + offset;
    }
    vector<bool> reconstruct() {
        vector<bool> res(N);
        queue<int> que;
        que.push(s);
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            if (s < N) res[s] = true;
            for (const auto &e : list[v]) {
                if (e.cap && !res[e.to]) {
                    res[e.to] = true;
                    que.push(e.to);
                }
            }
        }
        return res;
    }
    
private:
    // add edge
    Edge &get_rev_edge(const Edge &e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    Edge &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const Edge &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<Edge> get_edges() const {
        vector<Edge> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    void add_edge(int from, int to, T cap) {
        if (!cap) return;
        pos.emplace_back(from, (int)list[from].size());
        list[from].push_back(Edge((int)list[to].size(), from, to, cap));
        list[to].push_back(Edge((int)list[from].size() - 1, to, from, 0));
    }
    
    // Dinic's algorithm
    T dinic(int s, int t, T limit_flow) {
        T current_flow = 0;
        vector<int> level((int)list.size(), -1), iter((int)list.size(), 0);
        
        // Dinic BFS
        auto bfs = [&]() -> void {
            level.assign((int)list.size(), -1);
            level[s] = 0;
            queue<int> que;
            que.push(s);
            while (!que.empty()) {
                int v = que.front();
                que.pop();
                for (const Edge &e : list[v]) {
                    if (level[e.to] < 0 && e.cap > 0) {
                        level[e.to] = level[v] + 1;
                        if (e.to == t) return;
                        que.push(e.to);
                    }
                }
            }
        };
        
        // Dinic DFS
        auto dfs = [&](auto self, int v, T up_flow) {
            if (v == t) return up_flow;
            T res_flow = 0;
            for (int &i = iter[v]; i < (int)list[v].size(); ++i) {
                Edge &e = list[v][i], &re = get_rev_edge(e);
                if (level[v] >= level[e.to] || e.cap == 0) continue;
                T flow = self(self, e.to, min(up_flow - res_flow, e.cap));
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
                T flow = dfs(dfs, s, limit_flow - current_flow);
                if (!flow) break;
                current_flow += flow;
            }
        }
        return current_flow;
    };
    T dinic(int s, int t) {
        return dinic(s, t, numeric_limits<T>::max());
    }
};

    

/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void ARC_085_E() {
    int N;
    cin >> N;
    vector<long long> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    
    // i 個目の宝石を割らない: F, i 個目の宝石を割る: T とする
    TwoVariableSubmodularOpt<long long> tvs(N);
    for (int i = 0; i < N; ++i) {
        tvs.add_single_cost(i, -a[i], 0);
    }
    const long long INF = 1LL<<55;
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
            if ((j+1) % (i+1) == 0) {
                // i: T, j: F は禁止
                tvs.add_penalty(i, j, INF);
            }
        }
    }
    cout << -tvs.solve() << endl;
}


int main() {
    ARC_085_E();
}


