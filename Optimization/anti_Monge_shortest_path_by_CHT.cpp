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


// Dynamic Li Chao Tree for Pairwise Single Crossing Functions
/*
    Func: Pairwise Single Crossing 性を仮定
    注意点：i < j なる領域だけでなく、全領域での Pairwise Single Crossing 性が必要！！！

    - MIN: クエリ x の取りうる最小値
    - MAX: クエリ x の取りうる最大値 (f_i(MAX) がオーバーフローしないように注意）
    - INF: 最大値
 
    - insert (Func f_i): add f_i, O(log X)
    - query (x): min_i{ f_i(x) }, O(log X)
*/
template<class VAL> struct DynamicLiChaoTreeForPSC {
    using Func = function<VAL(VAL)>;
    struct IndexedFunc {
        Func f;
        int id;
        IndexedFunc(Func f, int id) : f(std::move(f)), id(id) {}
        VAL eval(VAL x) const { return f(x); };
    };
    struct Node {
        IndexedFunc func;
        unique_ptr<Node> left, right;
        explicit Node(IndexedFunc f) : func(std::move(f)) {}
    };
    
    const VAL MIN, MAX, INF = numeric_limits<VAL>::max() / 2;
    unique_ptr<Node> root;
    int next_id = 0;
    
    DynamicLiChaoTreeForPSC() {}
    DynamicLiChaoTreeForPSC(VAL MIN, VAL MAX, VAL INF = numeric_limits<VAL>::max() / 2) : MIN(MIN), MAX(MAX), INF(INF) {
        assert(MIN <= MAX);
    }

    // insert
    void insert(unique_ptr<Node> &p, VAL low, VAL high, IndexedFunc &f) {
        if (!p) {
            p = make_unique<Node>(std::move(f));
            return;
        }
        if (p->func.eval(low) <= f.eval(low) && p->func.eval(high) <= f.eval(high)) return;
        if (p->func.eval(low) >= f.eval(low) && p->func.eval(high) >= f.eval(high)) {
            p->func = std::move(f);
            return;
        }
        VAL mid = (low + high) / 2;
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
    pair<VAL, int> query(const unique_ptr<Node>& p, VAL low, VAL high, VAL x) const {
        if (!p) return {INF, -1};
        pair<VAL, int> res = {p->func.eval(x), p->func.id};
        if (low == high) return res;
        VAL mid = (low + high) / 2;
        if (x <= mid) return min(res, query(p->left, low, mid, x));
        else return min(res, query(p->right, mid + 1, high, x));
    }
    pair<VAL, int> query(VAL x) const {
        assert(MIN <= x && x <= MAX);
        return query(root, MIN, MAX, x);
    }
};

// find shortest path on DAG with Pairwise Single Crossing cost, by Li Chao Tree, O(N log N)
// vertex: 0, 1, 2, ..., N
// A(j, i) = dp[i] + f(i, j) must be Pairwise Single Crossing (must be in all 0 <= i, j <= N, including i >= j)
// Attention: Even if f is Pairwise Single Crossing, A is not necessarily Pairwise Single Crossing
// edges: (i, j), j >= start[i] (i + 1 <= start[i] <= N + 1), start is non-decreasing
template<class VAL, class FUNC> vector<pair<VAL, int>> PSCShortestPathByCHT(int N, const FUNC &f, const vector<int> &start) {
    assert((int)start.size() == N);
    for (int i = 0; i < N; i++) {
        assert(i + 1 <= start[i] && start[i] <= N + 1);
        if (i > 0) assert(start[i - 1] <= start[i]);
    }
    VAL INF = numeric_limits<VAL>::max() / 2;
    vector<pair<VAL, int>> dp(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
    DynamicLiChaoTreeForPSC<VAL> cht(1, N, INF);
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
template<class VAL, class FUNC> vector<pair<VAL, int>> PSCShortestPathByCHT(int N, const FUNC &f) {
    vector<int> start(N);
    for (int i = 0; i < N; i++) start[i] = i + 1;
    return PSCShortestPathByCHT<VAL>(N, f, start);
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
    auto res = PSCShortestPathByCHT<long long>(N, func, start);
    cout << res.back().first << endl;
}


int main() {
    ABC_480_G_CHT();
}