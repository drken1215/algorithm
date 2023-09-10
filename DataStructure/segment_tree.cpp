//
// セグメント木
//
// verified:
//   AtCoder Library Practice Contest J - Segment Tree (for max_right)
//     https://atcoder.jp/contests/practice2/tasks/practice2_j
//
//   ARC 008 D - タコヤキオイシクナール
//     https://beta.atcoder.jp/contests/arc008/tasks/arc008_4
//
//   ABC 281 E - Least Elements (for max_right)
//     https://atcoder.jp/contests/abc281/tasks/abc281_e
//
//   ABC 307 F - Virus 2 (for max_right)
//     https://atcoder.jp/contests/abc307/tasks/abc307_f
//


#include<bits/stdc++.h>
using namespace std;


// Segment Tree
template<class Monoid> struct SegTree {
    using Func = function<Monoid(Monoid, Monoid)>;

    // core member
    int N;
    Func OP;
    Monoid IDENTITY;
    
    // inner data
    int offset;
    vector<Monoid> dat;

    // constructor
    SegTree() {}
    SegTree(int n, const Func &op, const Monoid &identity) {
        init(n, op, identity);
    }
    SegTree(const vector<Monoid> &v, const Func &op, const Monoid &identity) {
        init((int)v.size(), op, identity);
        build(v);
    }
    void init(int n, const Func &op, const Monoid &identity) {
        N = n;
        OP = op;
        IDENTITY = identity;
        offset = 1;
        while (offset < N) offset *= 2;
        dat.assign(offset * 2, IDENTITY);
    }
    void init(const vector<Monoid> &v, const Func &op, const Monoid &identity) {
        init((int)v.size(), op, identity);
        build(v);
    }
    void build(const vector<Monoid> &v) {
        assert(N == (int)v.size());
        for (int i = 0; i < N; ++i) dat[i + offset] = v[i];
        for (int k = offset - 1; k > 0; --k) dat[k] = OP(dat[k*2], dat[k*2+1]);
    }
    int size() const {
        return N;
    }
    Monoid operator [] (int a) const { return dat[a + offset]; }
    
    // update A[a], a is 0-indexed, O(log N)
    void set(int a, const Monoid &v) {
        int k = a + offset;
        dat[k] = v;
        while (k >>= 1) dat[k] = OP(dat[k*2], dat[k*2+1]);
    }
    
    // get [a, b), a and b are 0-indexed, O(log N)
    Monoid prod(int a, int b) {
        Monoid vleft = IDENTITY, vright = IDENTITY;
        for (int left = a + offset, right = b + offset; left < right;
        left >>= 1, right >>= 1) {
            if (left & 1) vleft = OP(vleft, dat[left++]);
            if (right & 1) vright = OP(dat[--right], vright);
        }
        return OP(vleft, vright);
    }
    Monoid all_prod() { return dat[1]; }
    
    // get max r that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int max_right(const function<bool(Monoid)> f, int l = 0) {
        if (l == N) return N;
        l += offset;
        Monoid sum = IDENTITY;
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(OP(sum, dat[l]))) {
                while (l < offset) {
                    l = l * 2;
                    if (f(OP(sum, dat[l]))) {
                        sum = OP(sum, dat[l]);
                        ++l;
                    }
                }
                return l - offset;
            }
            sum = OP(sum, dat[l]);
            ++l;
        } while ((l & -l) != l);  // stop if l = 2^e
        return N;
    }

    // get min l that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int min_left(const function<bool(Monoid)> f, int r = -1) {
        if (r == 0) return 0;
        if (r == -1) r = N;
        r += offset;
        Monoid sum = IDENTITY;
        do {
            --r;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(OP(dat[r], sum))) {
                while (r < offset) {
                    r = r * 2 + 1;
                    if (f(OP(dat[r], sum))) {
                        sum = OP(dat[r], sum);
                        --r;
                    }
                }
                return r + 1 - offset;
            }
            sum = OP(dat[r], sum);
        } while ((r & -r) != r);
        return 0;
    }
    
    // debug
    friend ostream& operator << (ostream &s, const SegTree &seg) {
        for (int i = 0; i < seg.size(); ++i) {
            s << seg[i];
            if (i != seg.size()-1) s << " ";
        }
        return s;
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

// ACL practice J - Segment Tree
void ACL_practice_J() {
    int N, Q;
    cin >> N >> Q;
    vector<int> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    
    // セグ木の準備 (型: int, 演算方法: op, 単位元: -INF)
    const int INF = 1<<30;
    SegTree<int> seg(A, [&](int a, int b){ return max(a, b); }, -INF);
    
    // クエリ処理
    while (Q--) {
        int t;
        cin >> t;
        if (t == 1) {
            int X, V;
            cin >> X >> V;
            --X;
            
            // A[x] を V に update する処理
            seg.set(X, V);
        } else if (t == 2) {
            int L, R;
            cin >> L >> R;
            --L;
            
            // A の区間 [L, R) の最大値を求める処理
            cout << seg.prod(L, R) << endl;
        } else {
            int L, V;
            cin >> L >> V;
            --L;
            
            // x = seg.prod(L, r) として、f(x) = True となる最大の r を求める
            int res = seg.max_right([&](int x) -> bool { return V > x; }, L);
            
            // 求めたいのは V <= prod(L, x) となる最小の x
            cout << res + 1 << endl;
        }
    }
}

// ARC 008 D - タコヤキオイシクナール
void ARC_008_D() {
    // 入力
    long long N;
    int M;
    cin >> N >> M;
    vector<long long> p(M), pls;
    vector<double> a(M), b(M);
    for (int i = 0; i < M; ++i) {
        cin >> p[i] >> a[i] >> b[i];
        --p[i];
        pls.push_back(p[i]);
    }
    sort(pls.begin(), pls.end());
    pls.erase(unique(pls.begin(), pls.end()), pls.end());
    
    // セグメント木
    int NN = (int)pls.size();
    auto func = [&](pair<double,double> a, pair<double,double> b) {
        return make_pair(a.first * b.first, a.second * b.first + b.second);
    };
    pair<double, double> identity = make_pair(1, 0);
    SegTree<pair<double, double>> seg(NN, func, identity);

    // 処理
    double Min = 1.0, Max = 1.0;
    for (int i = 0; i < M; ++i) {
        int idx = lower_bound(pls.begin(), pls.end(), p[i]) - pls.begin();
        
        seg.set(idx, make_pair(a[i], b[i]));
        pair<double,double> res = seg.prod(0, NN);
        
        Min = min(Min, res.first + res.second);
        Max = max(Max, res.first + res.second);
    }
    cout << fixed << setprecision(10) << Min << endl << Max << endl;
}

// ABC 281 E
void ABC_281_E() {
    // 入力
    long long N, M, K;
    cin >> N >> M >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];

    // 座標圧縮を準備
    vector<pair<long long,int>> comp(N);
    for (int i = 0; i < N; ++i) comp[i] = {A[i], i};
    sort(comp.begin(), comp.end());
    vector<int> order(N);
    for (int i = 0; i < N; ++i) order[comp[i].second] = i;

    // セグメント木
    using Monoid = pair<int, long long>;
    auto func = [&](Monoid a, Monoid b) -> Monoid {
        return Monoid(a.first + b.first, a.second + b.second);
    };
    Monoid zero = Monoid(0, 0);
    SegTree<Monoid> seg(N, func, zero);

    // push と pop
    auto push = [&](long long x, int id) -> void {
        seg.set(id, Monoid(1, x));
    };
    auto pop = [&](long long x, int id) -> void {
        seg.set(id, zero);
    };
    auto get = [&]() -> long long {
        int r = seg.max_right([&](Monoid r) {
            return r.first <= K;
        });
        return seg.prod(0, r).second;
    };

    for (int i = 0; i < M; ++i) push(A[i], order[i]);
    for (int i = 0; i < N - M + 1; ++i) {
        cout << get() << " ";
        if (i+M < N) {
            push(A[i+M], order[i+M]);
            pop(A[i], order[i]);
        }
    }
    cout << endl;
}

// ABC 307 F
void ABC_307_F() {
    const long long INF = 1LL<<60;
    
    // 入力
    using Edge = pair<int, long long>;
    using Graph = vector<vector<Edge>>;
    int N, M, K, D;
    cin >> N >> M;
    Graph G(N);
    for (int i = 0; i < M; ++i) {
        int u, v, w;
        cin >> u >> v >> w;
        --u, --v;
        G[u].emplace_back(v, w);
        G[v].emplace_back(u, w);
    }
    cin >> K;
    vector<int> A(K);
    for (int i = 0; i < K; ++i) {
        cin >> A[i];
        --A[i];
    }
    cin >> D;
    vector<long long> X(D);
    for (int i = 0; i < D; ++i) cin >> X[i];
    
    // day >= start かつ X[day] >= x となる最小の day を求める
    SegTree<long long> seg(D, [&](long long a, long long b){return max(a,b);}, -INF);
    seg.build(X);
    auto first_okay_day = [&](long long x, int start) -> int {
        return seg.max_right([&](long long val) {return val < x;}, start);
    };

    // Dijkstra
    using Weight = pair<long long, long long>;
    using Node = pair<Weight, int>;
    priority_queue<Node, vector<Node>, greater<Node>> que;
    vector<Weight> dp(N, Weight(INF,INF));
    for (int i = 0; i < K; ++i) {
        que.push(Node({0,0}, A[i]));
        dp[A[i]] = {0,0};
    }
    while (!que.empty()) {
        auto [cur, v] = que.top();
        que.pop();
        
        if (cur > dp[v] || cur.first >= D) continue;
        for (auto e : G[v]) {
            Weight nex(cur.first, cur.second + e.second);
            if (cur.second + e.second > X[cur.first]) {
                nex = {first_okay_day(e.second, cur.first+1), e.second};
            }
            if (dp[e.first] > nex) {
                dp[e.first] = nex;
                que.push(Node(dp[e.first], e.first));
            }
        }
    }
    for (int v = 0; v < N; ++v) {
        if (dp[v].second > 0) ++dp[v].first;
        cout << (dp[v].first <= D ? dp[v].first : -1) << endl;
    }
}
    

int main() {
    ACL_practice_J();
    //ARC_008_D();
    //ABC_281_E();
    //ABC_307_F();
}

