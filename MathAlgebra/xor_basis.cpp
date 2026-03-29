//
// XOR 基底
//
// verified:
//   ZONeエナジー プログラミングコンテスト “HELLO SPACE” F - 出会いと別れ
//     https://atcoder.jp/contests/zone2021/tasks/zone2021_f
//


#include <bits/stdc++.h>
using namespace std;


#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl

// input
template<class T> istream& operator >> (istream &is, vector<T> &P)
{ for (int i = 0; i < P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, deque<T> &P)
{ for (int i = 0; i < P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, vector<vector<T>> &P)
{ for (int i = 0; i < P.size(); ++i) cin >> P[i]; return is; }

// output
template<class S, class T> ostream& operator << (ostream &s, const pair<S, T> &P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 2> &P)
{ return s << '<' << P[0] << "," << P[1] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 3> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 4> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << "," << P[3] << '>'; }
template<class T> ostream& operator << (ostream &s, const vector<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const deque<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const vector<vector<T>> &P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, const set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const multiset<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const unordered_set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const unordered_map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
void yes(bool a) { cout << (a ? "yes" : "no") << endl; }
void YES(bool a) { cout << (a ? "YES" : "NO") << endl; }
void Yes(bool a) { cout << (a ? "Yes" : "No") << endl; }
const vector<int> DX = {1, 0, -1, 0, 1, -1, 1, -1};
const vector<int> DY = {0, 1, 0, -1, 1, -1, -1, 1};



// xor basis
template<class T> vector<T> xor_basis(const vector<T> &vec) {
    vector<T> base, elim;
    for (auto v : vec) {
        T v2 = v;
        for (auto e : elim) v2 = min(v2, v2 ^ e);
        if (v2) base.emplace_back(v), elim.emplace_back(v2);
    }
    return base;
}


//------------------------------//
// Examples
//------------------------------//

// ZONeエナジー プログラミングコンテスト “HELLO SPACE” F - 出会いと別れ
struct UnionFind {
    // core member
    vector<int> par, nex;

    // constructor
    UnionFind() { }
    UnionFind(int N) : par(N, -1), nex(N) {
        init(N);
    }
    void init(int N) {
        par.assign(N, -1);
        nex.resize(N);
        for (int i = 0; i < N; ++i) nex[i] = i;
    }
    
    // core methods
    int root(int x) {
        if (par[x] < 0) return x;
        else return par[x] = root(par[x]);
    }
    
    bool same(int x, int y) {
        return root(x) == root(y);
    }
    
    bool merge(int x, int y, bool merge_technique = true) {
        x = root(x), y = root(y);
        if (x == y) return false;
        if (merge_technique) if (par[x] > par[y]) swap(x, y); // merge technique
        par[x] += par[y];
        par[y] = x;
        swap(nex[x], nex[y]);
        return true;
    }
    
    int size(int x) {
        return -par[root(x)];
    }
    
    // get group
    vector<int> group(int x) {
        vector<int> res({x});
        while (nex[res.back()] != x) res.push_back(nex[res.back()]);
        return res;
    }
    vector<vector<int>> groups() {
        vector<vector<int>> member(par.size());
        for (int v = 0; v < (int)par.size(); ++v) {
            member[root(v)].push_back(v);
        }
        vector<vector<int>> res;
        for (int v = 0; v < (int)par.size(); ++v) {
            if (!member[v].empty()) res.push_back(member[v]);
        }
        return res;
    }
    
    // debug
    friend ostream& operator << (ostream &s, UnionFind uf) {
        const vector<vector<int>> &gs = uf.groups();
        for (const vector<int> &g : gs) {
            s << "group: ";
            for (int v : g) s << v << " ";
            s << endl;
        }
        return s;
    }
};
void ZONe_F() {
    long long N, M;
    cin >> N >> M;
    vector<long long> all, A(M), can(N, 1);
    for (int i = 0; i < M; i++) cin >> A[i], can[A[i]] = 0;
    for (int i = 0; i < N; i++) if (can[i]) all.emplace_back(i);
    const auto &base = xor_basis(all);
    UnionFind uf(N);
    vector<pair<long long, long long>> res;
    for (auto v : base) {
        for (int i = 0; i < N; i++) {
            int j = i ^ v;
            if (!uf.same(i, j)) {
                uf.merge(i, j);
                res.emplace_back(i, j);
            }
        }
    }
    if (res.size() == N-1) {
        for (auto [x, y] : res) cout << x << " " << y << endl;
    } else cout << -1 << endl;
}


int main() {
    ZONe_F();
}