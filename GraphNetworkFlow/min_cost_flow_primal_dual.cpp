//
// min-cost flow (primal-dual)
//
// verified
//   AOJ Course GRL_6_B Network Flow - Minimum Cost Flow
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_6_B&lang=jp
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

template<class FLOWTYPE, class COSTTYPE> COSTTYPE MinCostFlow(Graph<FLOWTYPE, COSTTYPE> &G, int s, int t, FLOWTYPE f) {
    int n = (int)G.size();
    vector<COSTTYPE> pot(n, 0), dist(n, -1);
    vector<int> prevv(n), preve(n);
    COSTTYPE res = 0;
    while (f > 0) {
        priority_queue<pair<COSTTYPE,int>, vector<pair<COSTTYPE,int> >, greater<pair<COSTTYPE,int> > > que;
        dist.assign(n, -1);
        dist[s] = 0;
        que.push(make_pair(0,s));
        while(!que.empty()) {
            pair<COSTTYPE,int> p = que.top();
            que.pop();
            int v = p.second;
            if (dist[v] < p.first) continue;
            for (int i = 0; i < G[v].size(); ++i) {
                auto e = G[v][i];
                if (e.cap > 0 && (dist[e.to] < 0 || dist[e.to] > dist[v] + e.cost + pot[v] - pot[e.to])) {
                    dist[e.to] = dist[v] + e.cost + pot[v] - pot[e.to];
                    prevv[e.to] = v;
                    preve[e.to] = i;
                    que.push(make_pair(dist[e.to], e.to));
                }
            }
        }
        if (dist[t] < 0) return -1;
        for (int v = 0; v < n; ++v) pot[v] += dist[v];
        FLOWTYPE d = f;
        for (int v = t; v != s; v = prevv[v]) {
            d = min(d, G[prevv[v]][preve[v]].cap);
        }
        f -= d;
        res += pot[t] * d;
        for (int v = t; v != s; v = prevv[v]) {
            Edge<FLOWTYPE,COSTTYPE> &e = G[prevv[v]][preve[v]];
            Edge<FLOWTYPE,COSTTYPE> &re = G.redge(e);
            e.cap -= d;
            re.cap += d;
        }
    }
    return res;
}



int main() {
    int V, E, F; cin >> V >> E >> F;
    Graph<int,int> G(V);
    for (int i = 0; i < E; ++i) {
        int u, v, cap, cost; cin >> u >> v >> cap >> cost;
        G.addedge(u, v, cap, cost);
    }
    int s = 0, t = V-1;
    cout << MinCostFlow(G, s, t, F) << endl;
}  
