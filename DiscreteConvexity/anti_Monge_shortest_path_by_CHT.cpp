//
// anti-Monge 単一始点最短路問題 by Convex Hull Trick
//   頂点数 N+1 の DAG, 頂点 i, j 間のコスト f(i, j) が anti-Monge であることを仮定
//   O(N log N)
//
// verified
//   ABC 480 G - Zonal Score Minimization
//     https://atcoder.jp/contests/abc480/tasks/abc480_g
//


#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Template
//------------------------------//

using ll = long long;
using i128 = __int128_t;
using u128 = __uint128_t;
using pint = pair<int, int>;
using pll = pair<long long, long long>;
using tll = array<long long, 3>;
using fll = array<long long, 4>;
using vint = vector<int>;
using vll = vector<long long>;
using dint = deque<int>;
using dll = deque<long long>;
using vvint = vector<vector<int>>;
using vvll = vector<vector<long long>>;
using vpll = vector<pair<long long, long long>>;
template<class T> using min_priority_queue = priority_queue<T, vector<T>, greater<T>>;

template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }
template<class S, class T> inline auto maxll(S a, T b) { return max(ll(a), ll(b)); }
template<class S, class T> inline auto minll(S a, T b) { return min(ll(a), ll(b)); }
template<class T> auto max(const T &a) { return *max_element(a.begin(), a.end()); }
template<class T> auto min(const T &a) { return *min_element(a.begin(), a.end()); }
template<class T> auto argmax(const T &a) { return max_element(a.begin(), a.end()) - a.begin(); }
template<class T> auto argmin(const T &a) { return min_element(a.begin(), a.end()) - a.begin(); }
template<class T> auto accum(const vector<T> &a) { return accumulate(a.begin(), a.end(), T()); }
template<class T> auto accum(const deque<T> &a) { return accumulate(a.begin(), a.end(), T()); }

#define REP(i, a) for (long long i = 0; i < (long long)(a); i++)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); i++)
#define RREP(i, a) for (long long i = (a)-1; i >= (long long)(0); --i)
#define RREP2(i, a, b) for (long long i = (b)-1; i >= (long long)(a); --i)
#define EB emplace_back
#define PF push_front
#define PB push_back
#define MP make_pair
#define FI first
#define SE second
#define ALL(x) x.begin(), x.end()

// input stream
template<class T> istream& operator >> (istream &is, vector<T> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, deque<T> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, vector<vector<T>> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }

// output stream
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class S, class T> ostream& operator << (ostream &s, const pair<S, T> &P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 2> &P)
{ return s << '<' << P[0] << "," << P[1] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 3> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 4> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << "," << P[3] << '>'; }
template<class T> ostream& operator << (ostream &s, const vector<string> &P)
{ for (int i = 0; i < P.size(); ++i) { s << P[i] << endl; } return s; }
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


//------------------------------//
// Monge
//------------------------------//

// Convex Hull Trick
/*
    Func: Monge 性を仮定
    注意点：SMAWK や LARSCH とは異なり、i < j なる領域だけでなく、全領域での Monge 性が必要！！！

    - MIN: クエリ x の取りうる最小値
    - MAX: クエリ x の取りうる最大値 (f_i(MAX) がオーバーフローしないように注意）
    - INF: 最大値
 
    - insert (Func f_i): add f_i, O(log X)
    - query (x): min_i{ f_i(x) }, O(log X)
*/
template<class T> struct MongeCHT {
    using Func = function<T(T)>;
    struct IndexedFunc {
        Func f;
        int id;
        IndexedFunc(Func f, int id) : f(std::move(f)), id(id) {}
        T eval(T x) const { return f(x); };
    };
    struct Node {
        IndexedFunc func;
        unique_ptr<Node> left, right;
        explicit Node(IndexedFunc f) : func(std::move(f)) {}
    };
    
    const T MIN, MAX, INF = numeric_limits<T>::max() / 2;
    unique_ptr<Node> root;
    int next_id = 0;
    
    MongeCHT() {}
    MongeCHT(T MIN, T MAX, T INF = numeric_limits<T>::max() / 2) : MIN(MIN), MAX(MAX), INF(INF) {
        assert(MIN <= MAX);
    }

    // insert
    void insert(unique_ptr<Node> &p, T low, T high, IndexedFunc &f) {
        if (!p) {
            p = make_unique<Node>(std::move(f));
            return;
        }
        if (p->func.eval(low) <= f.eval(low) && p->func.eval(high) <= f.eval(high)) return;
        if (p->func.eval(low) >= f.eval(low) && p->func.eval(high) >= f.eval(high)) {
            p->func = std::move(f);
            return;
        }
        T mid = (low + high) / 2;
        if (p->func.eval(mid) > f.eval(mid)) swap(p->func, f);
        if (low == high) return;
        if (p->func.eval(low) >= f.eval(low)) insert(p->left, low, mid, f);
        else insert(p->right, mid + 1, high, f);
    }
    int insert(Func f) {
        int id = next_id++;
        IndexedFunc func(std::move(f), id);
        insert(root, MIN, MAX, func);
        return id;
    }
    void clear() {
        root.reset();
        next_id = 0;
    }

    // query
    pair<T, int> query(const unique_ptr<Node>& p, T low, T high, T x) const {
        if (!p) return {INF, -1};
        pair<T, int> res = {p->func.eval(x), p->func.id};
        if (low == high) return res;
        T mid = (low + high) / 2;
        if (x <= mid) return min(res, query(p->left, low, mid, x));
        else return min(res, query(p->right, mid + 1, high, x));
    }
    pair<T, int> query(T x) const {
        assert(MIN <= x && x <= MAX);
        return query(root, MIN, MAX, x);
    }
};

// find shortest path on DAG with Monge or anti-Monge cost, by Convex Hull Trick, O(N log N)
// vertex: 0, 1, 2, ..., N
// f(i, j) must be Monge or anti-Monge
// domain: j >= start[i], start is non-decreasing
template<class VAL, class FUNC> vector<pair<VAL, int>> MongeShortestPathByCHT(int N, const FUNC &f, const vector<int> &start) {
    assert((int)start.size() == N);
    VAL INF = numeric_limits<VAL>::max() / 2;
    vector<pair<VAL, int>> dp(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
    MongeCHT<VAL> cht(1, N, INF);
    dp[0].first = 0;
    int i = 0;
    for (int j = 1; j <= N; j++) {
        while (i < N && j >= start[i]) {
            auto f2 = [&, i](int j) -> long long { return dp[i].first + f(i, j); };
            cht.insert(f2);
            i++;
        }
        dp[j] = cht.query(j);
    };
    return dp;
}
template<class VAL, class FUNC> vector<pair<VAL, int>> MongeShortestPathByCHT(int N, const FUNC &f) {
    vector<int> start(N);
    for (int i = 0; i < N; i++) start[i] = i + 1;
    return MongeShortestPathByCHT<VAL>(N, f, start);
}


//------------------------------//
// Examples
//------------------------------//

// ABC 480 G - Zonal Score Minimization
template<class MeetSemiLattice> struct SparseTable {
    using Func = function<MeetSemiLattice(MeetSemiLattice, MeetSemiLattice)>;

    // core member
    Func OP = [](const MeetSemiLattice &l, const MeetSemiLattice &r) {
        return min(l, r);
    };
    vector<vector<MeetSemiLattice>> dat;
    vector<int> height;
    
    SparseTable() {}
    SparseTable(const vector<MeetSemiLattice> &vec) {
        init(vec);
    }
    SparseTable(const vector<MeetSemiLattice> &vec, const Func &op)  {
        init(vec, op);
    }
    void init(const vector<MeetSemiLattice> &vec) {
        int n = (int)vec.size(), h = 1;
        while ((1<<h) <= n) ++h;
        dat.assign(h, vector<MeetSemiLattice>(1<<h));
        height.assign(n+1, 0);
        for (int i = 2; i <= n; i++) height[i] = height[i>>1]+1;
        for (int i = 0; i < n; ++i) dat[0][i] = vec[i];
        for (int i = 1; i < h; ++i) {
            for (int j = 0; j < n; ++j)
                dat[i][j] = OP(dat[i-1][j], dat[i-1][min(j+(1<<(i-1)),n-1)]);
        }
    }
    void init(const vector<MeetSemiLattice> &vec, const Func &op) {
        OP = op;
        init(vec);
    }
    
    MeetSemiLattice get(int a, int b) {
        return OP(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};
void ABC_480_G_CHT() {
    int N;
    cin >> N;
    long long INF = numeric_limits<long long>::max() / 2;
    long long M = numeric_limits<long long>::max() / (N * 3);
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    SparseTable<long long> stmi(A);
    SparseTable<long long> stma(A, [&](long long p, long long q){return max(p, q);});
    auto func = [&](int l, int r) -> long long { 
        if (r < l + 2) return -M * ((l + 1 - r) * 2 + 1);
        else return stma.get(l, r) - stmi.get(l, r);
    };
    vector<int> start(N);
    for (int i = 0; i < N; i++) start[i] = i + 2;  // 辺 (i, i + 1) は禁止
    auto res = MongeShortestPathByCHT<long long>(N, func, start);
    cout << res.back().first << endl;
}


int main() {
    ABC_480_G_CHT();
}