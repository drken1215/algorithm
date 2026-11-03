//
// Dynamic Li Chao Tree
//   for Pairwise Single Crossing 性を満たす関数群 f_1, f_2, ..., f_N
//
// verified
//   Library Checker - Line Add Get Min
//     https://judge.yosupo.jp/problem/line_add_get_min
//
//   Library Checker - Segment Add Get Min
//     https://judge.yosupo.jp/problem/segment_add_get_min
//
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



// Dynamic Li Chao Tree for Pairwise Single Crossing Functions
/*
    Func: Pairwise Single Crossing 性を仮定
    注意点：i < j なる領域だけでなく、全領域での Pairwise Single Crossing 性が必要！！！

    x 座標のとりうる値の範囲: [XL, XR)
*/
template<class VAL, class FUNC, class X = long long> struct DynamicPSCLiChaoTree {
    using Func = function<VAL(X)>;
    struct IndexedFunc {
        Func f;
        int id;
        VAL eval(X x) const { return f(x); };
    };
    struct Node {
        optional<IndexedFunc> func;
        unique_ptr<Node> left, right;
    };
    
    // inner values
    X XL, XR;
    VAL INF = numeric_limits<VAL>::max() / 2;
    int next_id = 0;
    unique_ptr<Node> root;
    
    // constructors
    DynamicPSCLiChaoTree() {}
    DynamicPSCLiChaoTree(X xl, X xr, VAL inf = numeric_limits<VAL>::max() / 2) : XL(xl), XR(xr), INF(inf) {
        assert(xl < xr);
    }
    static bool better(const IndexedFunc &a, const IndexedFunc &b, X x) {
        return make_pair(a.eval(x), a.id) < make_pair(b.eval(x), b.id);
    }
    void clear() {
        root.reset();
        next_id = 0;
    }

    // insert(f), in O(log X)
    void insert(unique_ptr<Node> &p, X left, X right, IndexedFunc f) {
        if (!p) p = make_unique<Node>();
        if (!p->func) {
            p->func.emplace(f);
            return;
        }
        X mid = (left + right) / 2;
        bool left_better = better(f, *p->func, left);
        bool right_better = better(f, *p->func, right - 1);

        if (!left_better && !right_better) return;
        if (left_better && right_better) {
            p->func.emplace(f);
            return;
        }
        bool mid_better = better(f, *p->func, mid);
        if (mid_better) {
            IndexedFunc old = std::move(*p->func);
            p->func.emplace(f);
            if (left_better != mid_better) insert(p->left, left, mid, old);
            else insert(p->right, mid, right, old);
        } else {
            if (left_better != mid_better) insert(p->left, left, mid, f);
            else insert(p->right, mid, right, f);
        }
    }
    int insert(Func f) {
        int id = next_id++;
        insert(root, XL, XR, IndexedFunc{std::move(f), id});
        return id;
    }

    // insert_segment(l, r, f), in O(log^2 X)
    void insert_segment(unique_ptr<Node> &p, X left, X right, X l, X r, const IndexedFunc &f) {
        if (r <= left || right <= l) return;
        if (l <= left && right <= r) {
            insert(p, left, right, f);
            return;
        }
        if (!p) p = make_unique<Node>();
        X mid = (left + right) / 2;
        insert_segment(p->left, left, mid, l, r, f);
        insert_segment(p->right, mid, right, l, r, f);
    }
    int insert_segment(X l, X r, Func f) {
        assert(XL <= l && l <= r && r <= XR);
        int id = next_id++;
        if (l != r) insert_segment(root, XL, XR, l, r, IndexedFunc{std::move(f), id});
        return id;
    }

    // query(x), in O(log X)
    pair<VAL, int> query(X x) const {
        assert(XL <= x && x < XR);
        pair<VAL, int> res{INF, -1};
        const Node *p = root.get();
        X left = XL, right = XR;
        while (p) {
            if (p->func) {
                auto cur = make_pair(p->func->eval(x), p->func->id);
                if (res.second == -1 || cur < res) res = cur;
            }
            if (right - left == 1) break;
            X mid = (left + right) / 2;
            if (x < mid) {
                p = p->left.get();
                right = mid;
            } else {
                p = p->right.get();
                left = mid;
            }
        }
        return res;
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
    auto make_func = [&](int i) {
        return [&, i](long long j) -> VAL { return dp[i].first + f(i, j); };
    };
    DynamicPSCLiChaoTree<VAL, decltype(make_func(0))> cht(1, N + 1, INF);
    dp[0].first = 0;
    int i = 0;
    for (int j = 1; j <= N; j++) {
        while (i < N && j >= start[i]) {
            cht.insert(make_func(i));
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

// Library Checker - Line Add Get Min
void LibaryCheckerLineAddGetMin() {
    long long LX = -(1LL << 30), RX = 1LL << 30;
    auto make_func = [&](long long a, long long b) {
        return [&, a, b](long long x) -> long long {
            return a * x + b;
        };
    };
    DynamicPSCLiChaoTree<long long, decltype(make_func(0, 0))> cht(LX, RX);
    long long N, Q, typ, a, b, p;
    cin >> N >> Q;
    for (int i = 0; i < N; i++) {
        cin >> a >> b;
        cht.insert(make_func(a, b));
    }
    while (Q--) {
        cin >> typ;
        if (typ == 0) {
            cin >> a >> b;
            cht.insert(make_func(a, b));
        } else {
            cin >> p;
            auto [res, id] = cht.query(p);
            cout << res << endl;
        }
    }
}

// Library Checker - Segment Add Get Min
struct Line {
    long long a, b;
    long long operator () (long long x) {
        return a * x + b;
    }
};
void LibaryCheckerSegmentAddGetMin() {
    long long LX = -(1LL << 30), RX = 1LL << 30;
    auto make_func = [&](long long a, long long b) {
        return [&, a, b](long long x) -> long long {
            return a * x + b;
        };
    };
    DynamicPSCLiChaoTree<long long, decltype(make_func(0, 0))> cht(LX, RX);
    long long N, Q, typ, l, r, a, b, p;
    cin >> N >> Q;
    for (int i = 0; i < N; i++) {
        cin >> l >> r >> a >> b;
        cht.insert_segment(l, r, make_func(a, b));
    }
    while (Q--) {
        cin >> typ;
        if (typ == 0) {
            cin >> l >> r >> a >> b;
            cht.insert_segment(l, r, make_func(a, b));
        } else {
            cin >> p;
            auto [res, id] = cht.query(p);
            if (id == -1) cout << "INFINITY" << endl;
            else cout << res << endl;
        }
    }
}

// COLOCON 2018 Final C - スペースエクスプローラー高橋君
void COLOCON_2018_final_C() {
    long long N;
    cin >> N;
    vector<long long> a(N), res(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    auto make_func = [&](long long i) {
        return [&, i](long long x) -> long long {
            return (x - i) * (x - i) + a[i];
        };
    };
    DynamicPSCLiChaoTree<long long, decltype(make_func(0))> cht(0, 210000);
    for (long long i = 0; i < N; ++i) {
        cht.insert(make_func(i));
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
    //LibaryCheckerLineAddGetMin();
    //LibaryCheckerSegmentAddGetMin();
    //COLOCON_2018_final_C();
    EDPC_Z();
    //Codeforces_189_C();
    //yukicoder_705();
}