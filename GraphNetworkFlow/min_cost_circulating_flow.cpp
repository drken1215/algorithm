//
// Goldberg--Tarjan による cost-scaling を用いた最小費用循環流
//
// verified
//   AOJ Course GRL_6_B Network Flow - Minimum Cost Flow
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_6_B&lang=jp
//


#include <bits/stdc++.h>
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

template<class FLOWTYPE, class COSTTYPE> double newcost
(const Edge<FLOWTYPE, COSTTYPE> &e, const vector<double> &price) {
    return e.cost - price[e.from] + price[e.to];
}
    
template<class FLOWTYPE, class COSTTYPE> void DFSOfConstructGaux
(Graph<FLOWTYPE, COSTTYPE> &G, int v
, vector<bool> &visited, vector<double> &price) {
    visited[v] = true;
    for (int i = 0; i < G[v].size(); ++i) {
        Edge<FLOWTYPE, COSTTYPE> &e = G[v][i];
        if (e.cap > 0 && !visited[e.to] && newcost(e, price) < 0) {
            DFSOfConstructGaux(G, e.to, visited, price);
        }
    }
}

template<class FLOWTYPE, class COSTTYPE> void ConstructGaux
(Graph<FLOWTYPE, COSTTYPE> &G, double eps
, vector<FLOWTYPE> &balance, vector<double> &price) {
    vector<bool> visited(G.size(), false);
    for (int v = 0; v < G.size(); ++v) {
        if (balance[v] > 0)
            DFSOfConstructGaux(G, v, visited, price);
    }
    for (int v = 0; v < G.size(); ++v) {
        if (visited[v])
            price[v] += eps;
    }
}

template<class FLOWTYPE, class COSTTYPE> FLOWTYPE DFSOfAugmentBlockingFlow
(Graph<FLOWTYPE, COSTTYPE> &G, int v, FLOWTYPE flow
, vector<int> &iter, vector<FLOWTYPE> &balance, vector<double> &price) {
    if (balance[v] < 0) {
        FLOWTYPE dif = min(flow, -balance[v]);
        balance[v] += dif;
        return dif;
    }
    for (; iter[v] < G[v].size(); iter[v]++) {
        Edge<FLOWTYPE, COSTTYPE> &e = G[v][iter[v]];
        if (e.cap > 0 && newcost(e, price) < 0) {
            FLOWTYPE dif = DFSOfAugmentBlockingFlow(G, e.to, min(flow, e.cap), iter, balance, price);
            if (dif > 0) {
                e.cap -= dif;
                G.redge(e).cap += dif;
                return dif;
            }
        }
    }
    return 0;
}

template<class FLOWTYPE, class COSTTYPE> bool AugmentBlockingFlow
(Graph<FLOWTYPE, COSTTYPE> &G, vector<FLOWTYPE> &balance, vector<double> &price) {
    vector<int> iter(G.size(), 0);
    bool finish = true;
    for (int v = 0; v < G.size(); ++v) {
        FLOWTYPE flow;
        while (balance[v] > 0
               && (flow = DFSOfAugmentBlockingFlow(G, v, balance[v], iter, balance, price)) > 0)
            balance[v] -= flow;
        if (balance[v] > 0)
            finish = false;
    }
    if (finish) return true;
    else return false;
}

template<class FLOWTYPE, class COSTTYPE> COSTTYPE MinCostCirculation(Graph<FLOWTYPE, COSTTYPE> &G) {
    double eps = 0.0;
    vector<double> price(G.size(), 0.0);
    for (int v = 0; v < G.size(); ++v) {
        for (int i = 0; i < G[v].size(); ++i) {
            Edge<FLOWTYPE, COSTTYPE> &e = G[v][i];
            if (e.cap > 0) eps = max(eps, -(double)(e.cost));
        }
        price[v] = 0.0;
    }
    vector<FLOWTYPE> balance(G.size(), 0);
    while (eps * G.size() > 1) {
        eps /= 2;
        for (int v = 0; v < G.size(); ++v) {
            for (int i = 0; i < G[v].size(); ++i) {
                Edge<FLOWTYPE, COSTTYPE> &e = G[v][i];
                if (e.cap > 0 && newcost(e, price) < 0) {
                    balance[e.from] -= e.cap;
                    balance[e.to] += e.cap;
                    G.redge(e).cap += e.cap;
                    e.cap = 0;
                }
            }
        }
        while (true) {
            ConstructGaux(G, eps, balance, price);
            if (AugmentBlockingFlow(G, balance, price)) break;
        }
    }
    COSTTYPE res = 0;
    for (int v = 0; v < G.size(); ++v) {
        for (int i = 0; i < G[v].size(); ++i) {
            Edge<FLOWTYPE, COSTTYPE> &e = G[v][i];
            if (e.icap > e.cap)
                res += e.cost * (e.icap - e.cap);
        }
    }
    return res;
}



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// AOJ
void AOJ_Course_GRL_6_B() {
    const long long INF = 1<<29;  // 十分大きい値
    
    // 入力
    int V, E;
    long long F;
    cin >> V >> E >> F;
    Graph<long long, long long> G(V);
    for (int i = 0; i < E; ++i) {
        int u, v, cap, cost;
        cin >> u >> v >> cap >> cost;
        G.addedge(u, v, cap, cost);
    }
    
    // 流量 F を強制するために、t から s へコスト -INF の辺を張る
    int s = 0, t = V-1;
    G.addedge(t, s, F, -INF);
    
    // 最小費用循環流
    long long res = MinCostCirculation(G) + F * INF;
    
    // 流量 F を流せない場合
    if (res >= INF) cout << -1 << endl;
    else cout << res << endl;
}


int main() {
    AOJ_Course_GRL_6_B();
}


