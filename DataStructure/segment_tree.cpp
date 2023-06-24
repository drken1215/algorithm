//
// segment tree
//
// verified:
//   ARC 008 D - タコヤキオイシクナール
//     https://beta.atcoder.jp/contests/arc008/tasks/arc008_4
//
//   ABC 281 E - Least Elements
//     https://atcoder.jp/contests/abc281/tasks/abc281_e
//
//   ABC 307 F - Virus 2
//     https://atcoder.jp/contests/abc307/tasks/abc307_f
//


#include <bits/stdc++.h>
using namespace std;


// Segment Tree
template<class Monoid> struct SegTree {
    using Func = function<Monoid(Monoid, Monoid)>;

    // core member
    int SIZE;
    Func F;
    Monoid IDENTITY;
    
    // data
    int offset;
    vector<Monoid> dat;

    // constructor
    SegTree() {}
    SegTree(int n, const Func &f, const Monoid &identity)
    : SIZE(n), F(f), IDENTITY(identity) {
        offset = 1;
        while (offset < n) offset *= 2;
        dat.assign(offset * 2, IDENTITY);
    }
    void init(int n, const Func &f, const Monoid &identity) {
        SIZE = n;
        F = f;
        IDENTITY = identity;
        offset = 1;
        while (offset < n) offset *= 2;
        dat.assign(offset * 2, IDENTITY);
    }
    int size() const { return SIZE; }
    
    // set, a is 0-indexed //
    // build(): O(N)
    void set(int a, const Monoid &v) { dat[a + offset] = v; }
    void build() {
        for (int k = offset - 1; k > 0; --k)
            dat[k] = F(dat[k*2], dat[k*2+1]);
    }
    void build(const vector<Monoid> &vec) {
        for (int a = 0; a < vec.size() && a + offset < dat.size(); ++a)
            set(a, vec[a]);
        build();
    }
    
    // update a, a is 0-indexed, O(log N)
    void update(int a, const Monoid &v) {
        int k = a + offset;
        dat[k] = v;
        while (k >>= 1) dat[k] = F(dat[k*2], dat[k*2+1]);
    }
    
    // get [a, b), a and b are 0-indexed, O(log N)
    Monoid get(int a, int b) {
        Monoid vleft = IDENTITY, vright = IDENTITY;
        for (int left = a + offset, right = b + offset; left < right;
        left >>= 1, right >>= 1) {
            if (left & 1) vleft = F(vleft, dat[left++]);
            if (right & 1) vright = F(dat[--right], vright);
        }
        return F(vleft, vright);
    }
    Monoid get_all() { return dat[1]; }
    Monoid operator [] (int a) const { return dat[a + offset]; }
    
    // get max r that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int max_right(const function<bool(Monoid)> f, int l = 0) {
        if (l == SIZE) return SIZE;
        l += offset;
        Monoid sum = IDENTITY;
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(F(sum, dat[l]))) {
                while (l < offset) {
                    l = l * 2;
                    if (f(F(sum, dat[l]))) {
                        sum = F(sum, dat[l]);
                        ++l;
                    }
                }
                return l - offset;
            }
            sum = F(sum, dat[l]);
            ++l;
        } while ((l & -l) != l);  // stop if l = 2^e
        return SIZE;
    }

    // get min l that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int min_left(const function<bool(Monoid)> f, int r = -1) {
        if (r == 0) return 0;
        if (r == -1) r = SIZE;
        r += offset;
        Monoid sum = IDENTITY;
        do {
            --r;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(F(dat[r], sum))) {
                while (r < offset) {
                    r = r * 2 + 1;
                    if (f(F(dat[r], sum))) {
                        sum = F(dat[r], sum);
                        --r;
                    }
                }
                return r + 1 - offset;
            }
            sum = F(dat[r], sum);
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

// template chmax, chmin
template<class T> inline bool chmax(T& a, T b) { if (a < b) { a = b; return 1; } return 0; }
template<class T> inline bool chmin(T& a, T b) { if (a > b) { a = b; return 1; } return 0; }

// template debug stream
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

// ARC 008 D - タコヤキオイシクナール
void ARC_008_D() {
    // 入力
    long long N;
    int M;
    cin >> N >> M;
    vector<long long> p(M), pls;
    vector<double> a(M), b(M);
    for (int i = 0; i < M; ++i) {
        cin >> p[i] >> a[i] >> b[i]; --p[i];
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
        seg.update(idx, make_pair(a[i], b[i]));
        pair<double,double> res = seg.get(0, NN);
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
        seg.update(id, Monoid(1, x));
    };
    auto pop = [&](long long x, int id) -> void {
        seg.update(id, zero);
    };
    auto get = [&]() -> long long {
        int r = seg.max_right([&](Monoid r) {
            return r.first <= K;
        });
        return seg.get(0, r).second;
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
            if (chmin(dp[e.first], nex)) {
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
    //ARC_008_D();
    ABC_281_E();
    //ABC_307_F();
}

