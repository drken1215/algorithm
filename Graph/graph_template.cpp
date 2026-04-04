//
// グラフを「隣接リスト」として表現する 
//

#include <bits/stdc++.h>
using namespace std;


// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
template<class T = long long> struct Graph {
    int V, E;
    bool record_edge_index = false;
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    vector<unordered_map<int, int>> id;  // id[v][w] := the index of node w in G[v]

    // constructors
    Graph(int n = 0, int m = 0, bool rei = false) 
        : V(n), E(m), record_edge_index(rei)
        , list(n), reversed_list(n), id(n) { }
    void init(int n = 0, int m = 0, bool rei = false) {
        V = n, E = m, record_edge_index = rei;
        list.assign(n, vector<Edge<T>>());
        reversed_list.assign(n, vector<Edge<T>>());
        if (record_edge_index) id.assign(n, unordered_map<int, int>());
    }
    Graph(const Graph&) = default;
    Graph& operator = (const Graph&) = default;

    // getters
    vector<Edge<T>> &operator [] (int i) { return list[i]; }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    const vector<Edge<T>> &get_rev_edges(int i) const { return reversed_list[i]; }
    constexpr size_t size() const { return list.size(); }
    constexpr void clear() { V = 0; list.clear(); }
    constexpr void resize(int n) { V = n; list.resize(n); }
    Edge<T> &get_edge(int u, int v) {
        assert(record_edge_index);
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u][v] >= 0 && id[u][v] < list[u].size());
        return list[u][id[u][v]];
    }
    const Edge<T> &get_edge(int u, int v) const {
        assert(record_edge_index);
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u].at(v) >= 0 && id[u].at(v) < list[u].size());
        return list[u][id[u].at(v)];
    }

    // add edge
    void add_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size(); 
        list[from].push_back(Edge(from, to, val));
        reversed_list[to].push_back(Edge(to, from, val));
    }
    void add_bidirected_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size(); 
        list[from].push_back(Edge(from, to, val));
        reversed_list[from].push_back(Edge(from, to, val));
        if (from != to) {
            if (record_edge_index) id[to][from] = (int)list[to].size(); 
            list[to].push_back(Edge(to, from, val));
            reversed_list[to].push_back(Edge(to, from, val));
        }
    }

    // input / output
    friend istream& operator >> (istream &is, Graph &G) {
        for (int i = 0; i < G.E; i++) {
            int u, v;
            is >> u >> v, u--, v--;
            G.add_bidirected_edge(u, v);
        }
        return is;
    }
    friend ostream &operator << (ostream &os, const Graph &G) {
        os << endl;
        for (int i = 0; i < (int)G.size(); ++i) {
            os << i << " -> ";
            for (int j = 0; j < (int)G[i].size(); j++) {
                if (j) os << ", ";
                os << G[i][j].to << "(" << G[i][j].val << ")";
            }
            os << endl;
        }
        return os;
    }
};


//------------------------------//
// Examples
//------------------------------//

int main() {

}