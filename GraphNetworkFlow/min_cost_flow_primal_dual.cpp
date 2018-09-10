//
// max-flow (Dinic's algorithm)
//
// verified
//   AOJ Course GRL_6_A Network Flow - Maximum Flow
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_6_A&lang=jp
//


#include <iostream>
#include <vector>
#include <queue>
using namespace std;

// edge class (for network-flow)
template<class FLOWTYPE, class COSTTYPE> struct Edge {
    int rev, from, to;
    FLOWTYPE cap, icap;
    COSTTYPE cost;
    Edge(int r, int f, int t, FLOWTYPE ca, COSTTYPE co) : rev(r), from(f), to(t), cap(ca), icap(ca), cost(co) {}
    friend ostream& operator << (ostream& s, const Edge& E) {
        if (E.cap > 0) return s << E.from << "->" << E.to << '(' << E.cap << ',' << E.cost << ')';
        else return s;
    }
};

// graph class (for network-flow)
template<class FLOWTYPE, class COSTTYPE> struct Graph {
    vector<vector<Edge<FLOWTYPE, COSTTYPE> > > list;
    
    Graph(int n = 0) : list(n) { }
    void init(int n = 0) { list.clear(); list.resize(n); }
    void reset() { for (int i = 0; i < (int)list.size(); ++i) for (int j = 0; j < list[i].size(); ++j) list[i][j].cap = list[i][j].icap; }
    inline vector<Edge<FLOWTYPE, COSTTYPE> >& operator [] (int i) { return list[i]; }
    inline const size_t size() const { return list.size(); }
    
    inline Edge<FLOWTYPE, COSTTYPE> &redge(const Edge<FLOWTYPE, COSTTYPE> &e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    
    void addedge(int from, int to, FLOWTYPE cap, COSTTYPE cost) {
        list[from].push_back(Edge<FLOWTYPE, COSTTYPE>((int)list[to].size(), from, to, cap, cost));
        list[to].push_back(Edge<FLOWTYPE, COSTTYPE>((int)list[from].size() - 1, to, from, 0, -cost));
    }
    
    void add_undirected_edge(int from, int to, FLOWTYPE cap, COSTTYPE cost) {
        list[from].push_back(Edge<FLOWTYPE, COSTTYPE>((int)list[to].size(), from, to, cap, cost));
        list[to].push_back(Edge<FLOWTYPE, COSTTYPE>((int)list[from].size() - 1, to, from, cap, cost));
    }

    /* 
    // debug
    friend ostream& operator << (ostream& s, const Graph& G) {
        s << endl; for (int i = 0; i < G.V; ++i) { s << i << " : " << G.list[i] << endl; }return s;
    }
    */
};

template<class FLOWTYPE, class COSTTYPE> struct MinCostFlow {
    const FLOWTYPE INF = 1<<30; // to be set
    vector<int> level, iter;

    Dinic() { }
    void dibfs(Graph<FLOWTYPE> &G, int s) {
        level.assign((int)G.size(), -1);
        level[s] = 0;
        queue<int> que;
        que.push(s);
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            for (int i = 0; i < G[v].size(); ++i) {
                Edge<FLOWTYPE> &e = G[v][i];
                if (level[e.to] < 0 && e.cap > 0) {
                    level[e.to] = level[v] + 1;
                    que.push(e.to);
                }
            }
        }
    }
    
    FLOWTYPE didfs(Graph<FLOWTYPE> &G, int v, int t, FLOWTYPE f) {
        if (v == t) return f;
        for (int &i = iter[v]; i < G[v].size(); ++i) {
            Edge<FLOWTYPE> &e = G[v][i], &re = G.redge(e);
            if (level[v] < level[e.to] && e.cap > 0) {
                FLOWTYPE d = didfs(G, e.to, t, min(f, e.cap));
                if (d > 0) {
                    e.cap -= d;
                    re.cap += d;
                    return d;
                }
            }
        }
        return 0;
    }
    
    FLOWTYPE solve(Graph<FLOWTYPE> &G, int s, int t) {
        level.assign((int)G.size(), -1); iter.assign((int)G.size(), 0);
        FLOWTYPE res = 0;
        while (true) {
            dibfs(G, s);
            if (level[t] < 0) return res;
            for (int i = 0; i < (int)iter.size(); ++i) iter[i] = 0;
            FLOWTYPE flow = 0;
            while ((flow = didfs(G, s, t, INF)) > 0) {
                res += flow;
            }
        }
    }
};



int main() {
    int V, E; cin >> V >> E;
    Graph<int> G(V);
    for (int i = 0; i < E; ++i) {
        int u, v, c; cin >> u >> v >> c;
        G.addedge(u, v, c);
    }
    Dinic<int> di;
    int s = 0, t = V-1;
    cout << di.solve(G, s, t) << endl;
}  
