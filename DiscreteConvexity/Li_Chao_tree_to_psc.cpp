//
// Dynamic Li Chao Tree
//   for Pairwise Single Crossing 性を満たす関数群 f_1, f_2, ..., f_N
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君 (for Monge CHT)
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//
//   AtCoder EDPC Z - Frog 3 (for simple LARSCH, LARSCH, monge CHT)
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//
//   Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry (for simple LARSCH, LARSCH, monge CHT)
//     https://codeforces.com/contest/319/problem/C 
// 
//   yukicoder No.705 ゴミ拾い Hard (for simple LARSCH, LARSCH, monge CHT)
//     https://yukicoder.me/problems/no/705 
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
template<class VAL> struct DynamicPSCLiChaoTree {
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
    
    DynamicPSCLiChaoTree() {}
    DynamicPSCLiChaoTree(VAL MIN, VAL MAX, VAL INF = numeric_limits<VAL>::max() / 2) : MIN(MIN), MAX(MAX), INF(INF) {
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
    DynamicPSCLiChaoTree<VAL> cht(1, N, INF);
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

// COLOCON 2018 Final C - スペースエクスプローラー高橋君
void COLOCON_2018_final_C() {
    long long N;
    cin >> N;
    vector<long long> a(N), res(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    DynamicLiChaoTreeForPSC<long long> cht(0, 210000, 1LL<<60);
    for (long long i = 0; i < N; ++i) {
        auto func = [i, &a](long long x) -> long long { return (x - i) * (x - i) + a[i]; };
        cht.insert(func);
    }
    for (long long i = 0; i < N; ++i) res[i] = cht.query(i).first;
    for (int i = 0; i < N; ++i) cout << res[i] << endl;
}

// AtCoder EDPC Z - Frog 3
/*
    H は単調増加数列
    chmin(dp[j], dp[i] + (H[j] - H[i])^2 + C)
    i -> j のコスト：(H[j] - H[i])^2 ...... 差の凸関数は Monge
    スタート: 0, ゴール: N-1
*/
void EDPC_Z() {
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N);
    for (long long i = 0; i < N; i++) cin >> H[i];
    auto func = [&](int i, int j) -> long long {
        return (H[j] - H[i]) * (H[j] - H[i]) + C;
    };
    const auto &res = PSCShortestPathByCHT<long long>(N-1, func);
    cout << res[N-1].first << endl;
}

// Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
/*
    A: 単調増加, B: 単調減少, ともに長さ N
    i -> j のコストが、B[i] × A[j] で与えられる ..... 単調増加 × 単調減少は Monge
    スタート: 0, ゴール: N-1
*/
void Codeforces_189_C() {
    long long N;
    cin >> N;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> B[i];
    auto func = [&](int i, int j) -> long long {
        return B[i] * A[j];
    };
    auto res = PSCShortestPathByCHT<long long>(N-1, func);
    cout << res[N-1].first << endl;
}

// yukicoder No.705 ゴミ拾い Hard
/*
    A, X, Y: N 個
    これらを区間に分割していく
    　dp[j] = min_{0 ≦ i < j}(dp[i] + |A[j-1] - X[i]|^3 + |-Y[i]|^3)
    i -> j のコスト：|A[j-1] - X[i]|^3 + |-Y[i]|^3 ...... 差の凸関数 (Monge) + 縞々 (Monge) -> Monge
    スタート: 0, ゴール: N
*/
void yukicoder_705() {
    int N;
    cin >> N;
    vector<long long> A(N), X(N), Y(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> X[i];
    for (int i = 0; i < N; i++) cin >> Y[i];
    auto func = [&](int i, int j) -> long long {
        long long dx = abs(A[j-1] - X[i]), dy = abs(Y[i]);
        return dx * dx * dx + dy * dy * dy;
    };
    auto res = PSCShortestPathByCHT<long long>(N, func);
    cout << res[N].first << endl;
}


int main() {
    //COLOCON_2018_final_C();
    EDPC_Z();
    //Codeforces_189_C();
    //yukicoder_705();
}