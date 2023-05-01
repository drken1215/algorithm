//
// 最小費用最大流 (正辺のみ)
//
// verified
//   AOJ 1088 School Excursion
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1088
//


#include <bits/stdc++.h>
using namespace std;
using pint = pair<int, int>;
using pll = pair<long long, long long>;
template<class T> inline bool chmax(T& a, T b) { if (a < b) { a = b; return 1; } return 0; }
template<class T> inline bool chmin(T& a, T b) { if (a > b) { a = b; return 1; } return 0; }

#define REP(i, n) for (long long i = 0; i < (long long)(n); ++i)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); ++i)
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, deque<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, vector<vector<T> > P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, set<T> P)
{ for(auto it : P) { s << "<" << it << "> "; } return s; }
template<class T1, class T2> ostream& operator << (ostream &s, map<T1,T2> P)
{ for(auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }


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

// min-cost max-flow (by primal-dual), maxf is the upper-limit of flow
template<class FLOWTYPE, class COSTTYPE> pair<FLOWTYPE, COSTTYPE> MinCostMaxFlow
 (Graph<FLOWTYPE, COSTTYPE> &G, int s, int t, FLOWTYPE maxf) {
    int n = (int)G.size();
    vector<COSTTYPE> pot(n, 0), dist(n, -1);
    vector<int> prevv(n), preve(n);
    FLOWTYPE f = maxf;
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
        if (dist[t] < 0) return {maxf - f, res};
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
    return {maxf, res};
}

int main() {
    int N, maxF;
    while (cin >> N, N) {
        // 入力
        vector<int> M(N-1);
        vector<vector<int>> X(N-1), Y(N-1), C(N-1);
        
        // (駅, 時刻) → 入って来る側の頂点番号
        map<pint, int> index;
        int vertex_num = 0;
        
        // 各駅の到着時刻と発車時刻
        vector<set<int>> intime(N), outtime(N);
        
        // 入力受け取り
        for (int i = 0; i < N-1; ++i) {
            cin >> M[i];
            X[i].resize(M[i]), Y[i].resize(M[i]), C[i].resize(M[i]);
            for (int j = 0; j < M[i]; ++j) {
                cin >> X[i][j] >> Y[i][j] >> C[i][j];
                outtime[i].insert(X[i][j]);
                intime[i+1].insert(Y[i][j]);
                if (!index.count(pint(i, X[i][j]))) {
                    index[pint(i, X[i][j])] = vertex_num;
                    vertex_num += 2;
                }
                if (!index.count(pint(i+1, Y[i][j]))) {
                    index[pint(i+1, Y[i][j])] = vertex_num;
                    vertex_num += 2;
                }
            }
        }
        cin >> maxF;
        
        // グラフを構築
        int source = vertex_num++, sink = vertex_num++;
        Graph<int, int> G(vertex_num);
        const int INF = maxF;  // 理論上流せる上限

        // ソース → 最初の駅
        for (auto t : outtime[0]) {
            G.addedge(source, index[pint(0, t)] + 1, INF, 0);
        }
        
        // 最後の駅 → シンク
        for (auto t : intime[N-1]) {
            G.addedge(index[pint(N-1, t)] + 1, sink, INF, 0);
        }
        
        // 各駅の intime (容量 1) を分裂
        for (int i = 0; i < N; ++i) {
            for (auto t : intime[i]) {
                G.addedge(index[pint(i, t)], index[pint(i, t)] + 1, 1, 0);
            }
        }
        
        // 各駅の intime -> outtime
        for (int i = 0; i < N; ++i) {
            for (auto t1 : intime[i]) {
                for (auto t2 : outtime[i]) {
                    if (t1 <= t2) {
                        G.addedge(index[pint(i, t1)] + 1, index[pint(i, t2)] + 1, INF, 0);
                    }
                }
            }
        }
        
        // 各電車
        for (int i = 0; i < N-1; ++i) {
            for (int j = 0; j < X[i].size(); ++j) {
                G.addedge(index[pint(i, X[i][j])] + 1, index[pint(i+1, Y[i][j])], 1, C[i][j]);
            }
        }

        // 最小費用最大流を流す
        auto res = MinCostMaxFlow(G, source, sink, maxF);
        cout << res.first << " " << res.second << endl;
    }
}

