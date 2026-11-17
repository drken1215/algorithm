#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// Edge Class
template<class VAL = long long> struct Edge {
    using value_type = VAL;

    // inner values
    int from, to;
    VAL val;

    // constructors
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, VAL v = 1) : from(f), to(t), val(v) {}
    Edge(const Edge&) = default;
    Edge& operator = (const Edge&) = default;

    // getters
    VAL get_val() const { return val; }
    VAL get_cost() const { return val; }
    VAL get_weight() const { return val; }

    // output
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
template<class VAL = long long> struct Graph {
    using value_type = VAL;
    using edge_type = Edge<VAL>;

    // inner values
    int V;
    bool record_edge_index = false;
    vector<vector<Edge<VAL>>> list;
    vector<unordered_map<int, int>> id;  // id[v][w] := the index of node w in G[v]

    // constructors
    Graph(int n = 0, bool rei = false) {
        init(n, rei);
    }
    void init(int n = 0, bool rei = false) {
        V = n;
        record_edge_index = rei;
        list.assign(n, vector<Edge<VAL>>());
        if (record_edge_index) id.assign(n, unordered_map<int, int>());
    }
    Graph(const Graph&) = default;
    Graph& operator = (const Graph&) = default;
    void resize(int n) { 
        V = n;
        list.resize(n);
    }
    void assign(int n) {
        V = n;
        list.assign(n, vector<Edge<VAL>>());
    }
    void clear() {
        V = 0;
        list.clear();
    }
    void clear(int v) {
        assert(v >= 0 && v < (int)list.size());
        list[v].clear();
    }

    // getters
    vector<Edge<VAL>> &operator [] (int i) { return list[i]; }
    const vector<Edge<VAL>> &operator [] (int i) const { return list[i]; }
    size_t size() const { return list.size(); }
    Edge<VAL> &get_edge(int u, int v) {
        assert(record_edge_index);
        assert(u >= 0 && u < (int)list.size() && v >= 0 && v < (int)list.size());
        assert(id[u].count(v) && id[u][v] >= 0 && id[u][v] < (int)list[u].size());
        return list[u][id[u][v]];
    }
    const Edge<VAL> &get_edge(int u, int v) const {
        assert(record_edge_index);
        assert(u >= 0 && u < (int)list.size() && v >= 0 && v < (int)list.size());
        assert(id[u].count(v) && id[u].at(v) >= 0 && id[u].at(v) < (int)list[u].size());
        return list[u][id[u].at(v)];
    }

    // add edge
    void add_edge(int from, int to, VAL val = 1) {
        assert(0 <= from && from < (int)list.size() && 0 <= to && to < (int)list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size(); 
        list[from].push_back(Edge(from, to, val));
    }
    void add_bidirected_edge(int from, int to, VAL val = 1) {
        assert(0 <= from && from < (int)list.size() && 0 <= to && to < (int)list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size();
        list[from].push_back(Edge(from, to, val));
        if (from != to) {
            if (record_edge_index) id[to][from] = (int)list[to].size(); 
            list[to].push_back(Edge(to, from, val));
        }
    }

    // input (only tree-case)
    friend istream& operator >> (istream &is, Graph &G) {
        for (int i = 0; i < G.V - 1; i++) {
            int u, v;
            is >> u >> v, u--, v--;
            G.add_bidirected_edge(u, v);
        }
        return is;
    }

    // output
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
