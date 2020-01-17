//
// min-cost flow (primal-dual, negative edges are ok but negative cycles are ng)
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
    int rev, from, to, id;
    FLOWTYPE cap, icap;
    COSTTYPE cost;
    Edge(int r, int f, int t, FLOWTYPE ca, COSTTYPE co, int id = -1) :
        rev(r), from(f), to(t), cap(ca), icap(ca), cost(co), id(id) {}
    friend ostream& operator << (ostream& s, const Edge& E) {
        if (E.cap > 0)
            return s << E.from << "->" << E.to <<
                '(' << E.cap << ',' << E.cost << ')';
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
    
    void addedge(int from, int to, FLOWTYPE cap, COSTTYPE cost, int id = -1) {
        list[from].push_back(Edge<FLOWTYPE, COSTTYPE>((int)list[to].size(), from, to, cap, cost, id));
        list[to].push_back(Edge<FLOWTYPE, COSTTYPE>((int)list[from].size() - 1, to, from, 0, -cost));
    }
    
    void add_undirected_edge(int from, int to, FLOWTYPE cap, COSTTYPE cost, int id = -1) {
        list[from].push_back(Edge<FLOWTYPE, COSTTYPE>((int)list[to].size(), from, to, cap, cost, id));
        list[to].push_back(Edge<FLOWTYPE, COSTTYPE>((int)list[from].size() - 1, to, from, cap, cost, id));
    }

    friend ostream& operator << (ostream& s, const Graph& G) {
        s << endl;
        for (int i = 0; i < G.size(); ++i) {
            s << i << ":";
            for (auto e : G.list[i]) s << " " << e;
            s << endl;
        }
        return s;
    }   
};

// min-cost flow (by primal-dual)
template<class FLOWTYPE, class COSTTYPE> COSTTYPE MinCostFlow(Graph<FLOWTYPE, COSTTYPE> &G, int s, int t, FLOWTYPE f) {
    int n = (int)G.size();
    vector<COSTTYPE> dist(n, -1);
    vector<int> prevv(n), preve(n), seen(n, 0);
    COSTTYPE res = 0;
    while (f > 0) {
        seen.assign(n, 0);
        dist[s] = 0;
        seen[s] = true;
        while (true) {
            bool update = false;
            for (int v = 0; v < n; ++v) {
                if (!seen[v]) continue;
                for (int i = 0; i < G[v].size(); ++i) {
                    Edge<FLOWTYPE, COSTTYPE> &e = G[v][i];
                    if (e.cap > 0 && (!seen[e.to] || dist[e.to] > dist[v] + e.cost)) {
                        dist[e.to] = dist[v] + e.cost;
                        prevv[e.to] = v;
                        preve[e.to] = i;
                        seen[e.to] = true;
                        update = true;
                    }
                }
            }
            if (!update) break;
        } 
        if (!seen[t]) return -1;
        FLOWTYPE d = f;
        for (int v = t; v != s; v = prevv[v]) {
            d = min(d, G[prevv[v]][preve[v]].cap);
        }
        f -= d;
        res += dist[t] * d;
        for (int v = t; v != s; v = prevv[v]) {
            Edge<FLOWTYPE, COSTTYPE> &e = G[prevv[v]][preve[v]];
            Edge<FLOWTYPE, COSTTYPE> &re = G.redge(e);
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
