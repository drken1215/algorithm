//
// 区間 (と値) を Set で管理する構造体
//
// {[1, 3), (7, 11)} に対して erase(5, 9) をするときの挙動:
//     this version: del(5, 7)  
//     another version: del(7, 11), add(5, 9) 
//
// verified
//   第六回 アルゴリズム実技検定 M - 等しい数 (for update)
//     https://atcoder.jp/contests/past202104-open/tasks/past202104_m
//
//   RUPC 2018 G - Elevator (for insert, same)
//     https://onlinejudge.u-aizu.ac.jp/problems/2880 
//
//   AtCoder ABC 255 Ex - Range Harvest Query (for update)
//     https://atcoder.jp/contests/abc255/tasks/abc255_h 
//
//   yukicoder No.674 n連勤 (for insert)
//     https://yukicoder.me/problems/no/674 
//
//   AtCoder ABC 330 E - Mex and Update (for insert, erase, mex)
//     https://atcoder.jp/contests/abc330/tasks/abc330_e 
//
//   第五回 アルゴリズム実技検定 N - 旅行会社 (for insert, erase, lower_bound)
//     https://atcoder.jp/contests/past202012-open/tasks/past202012_n
//
//   Code Festival 2015 予選 B D - マスと駒と色塗り (for insert, lower_bound)
//     https://atcoder.jp/contests/code-festival-2015-qualb/tasks/codefestival_2015_qualB_d
//
//   CPSCO 2019 Session 1 E - Exclusive OR Queries (for insert with del, lower_bound)
//     https://atcoder.jp/contests/cpsco2019-s1/tasks/cpsco2019_s1_e
//
//   yukicoder No.3017 交互浴 (for update)
//     https://yukicoder.me/problems/no/3017 
//


#include <bits/stdc++.h>
using namespace std;


// Interval Set
// T: type of range, VAL: data type
template<class T, class VAL = long long> struct IntervalSet {
    struct Node {
        T l, r;
        VAL val;
        Node(const T &l, const T &r, const VAL &val) : l(l), r(r), val(val) {}
        constexpr bool operator < (const Node &rhs) const {
            if (l != rhs.l) return l < rhs.l;
            else return r < rhs.r;
        }
        friend ostream& operator << (ostream &s, const Node &e) {
            return s << "([" << e.l << ", " << e.r << "): " << e.val << ")";
        }
    };

    // internal values
    const VAL identity;
    set<Node> S;

    // constructor
    IntervalSet(const VAL &identity = VAL()) : identity(identity) {}
    IntervalSet(const vector<VAL> &v, const VAL &identity = VAL()) : identity(identity) {
        vector<Node> vec;
        for (int l = 0; l < (int)v.size();) {
            int r = l;
            while (r < (int)v.size() && v[r] == v[l]) r++;
            vec.emplace_back(l, r, v[l]);
            l = r;
        }
        S = set<Node>(vec.begin(), vec.end());
    }

    // get the basic iterators
    constexpr typename set<Node>::iterator begin() { return S.begin(); }
    constexpr typename set<Node>::iterator end() { return S.end(); }

    // get the iterator of interval which contains p
    // not exist -> S.end()
    constexpr typename set<Node>::iterator get(const T &p) {
        auto it = S.upper_bound(Node(p, numeric_limits<T>::max(), 0));
        if (it == S.begin()) return S.end();
        it = prev(it);
        if (it->l <= p && p < it->r) return it;
        else return S.end();
    }

    // get the leftist iterator of interval which contains value >= p
    constexpr typename set<Node>::iterator lower_bound(const T &p) {
        auto it = get(p);
        if (it != S.end()) return it;
        return S.upper_bound(Node(p, numeric_limits<T>::max(), 0));
    }
    
    // exist the interval which contains p: true, [l, r): true
    constexpr bool covered(const T &p) {
        auto it = get(p);
        if (it != S.end()) return true;
        else return false;
    }
    constexpr bool covered(const T &l, const T &r) {
        assert(l <= r);
        if (l == r) return true;
        auto it = get(l);
        if (it != S.end() && r <= it->r) return true;
        else return false;
    }

    // is p, q in same interval?
    constexpr bool same(const T &p, const T &q) {
        if (!covered(p) || !covered(q)) return false;
        return get(p) == get(q);
    }

    // get the value of interval which contains p
    // not exist -> identity
    constexpr VAL get_val(const T &p) {
        auto it = get(p);
        if (it != S.end()) return it->val;
        else return identity;
    }
    VAL operator [] (const T &p) const {
        return get_val(p);
    }

    // get mex (>= p)
    constexpr T get_mex(const T &p = 0) {
        auto it = S.upper_bound(Node(p, numeric_limits<T>::max(), 0));
        if (it == S.begin()) return p;
        it = prev(it); 
        if (it->l <= p && p < it->r) return it->r;
        else return p;
    }

    // update [l, r) with value val / insert [l, r)
    // del: reflect effects of interval-delete
    // add: reflect effects of interval add
    template<class ADDFUNC, class DELFUNC> void update(T l, T r, const VAL &val, const ADDFUNC &add, const DELFUNC &del) {
        auto it = S.lower_bound(Node(l, 0, val));
        while (it != S.end() && it->l <= r) {
            if (it->l == r) {
                if (it->val ==val) {
                    del(r, it->r, val);
                    r = it->r;
                    it = S.erase(it);
                }
                break;
            }
            if (it->r <= r) {
                del(it->l, it->r, it->val);
                it = S.erase(it);
            } else {
                if (it->val == val) {
                    r = it->r;
                    del(it->l, it->r, it->val);
                    it = S.erase(it);
                } else {
                    del(it->l, r, it->val);
                    Node node = *it;
                    it = S.erase(it);
                    it = S.emplace_hint(it, r, node.r, node.val);
                }
            }
        }
        if (it != S.begin()) {
            it = prev(it);
            if (it->r == l) {
                if (it->val == val) {
                    del(it->l, it->r, it->val);
                    l = it->l;
                    it = S.erase(it);
                }
            } else if (l < it->r) {
                if (it->val == val) {
                    del(it->l, it->r, it->val);
                    l = min(l, it->l);
                    r = max(r, it->r);
                    it = S.erase(it);
                } else {
                    if (r < it->r) {
                        it = S.emplace_hint(next(it), r, it->r, it->val);
                        it = prev(it);
                    }
                    del(l, min(r, it->r), it->val);
                    Node node = *it;
                    it = S.erase(it);
                    it = S.emplace_hint(it, node.l, l, node.val);
                }
            }
        }
        if (it != S.end()) it = next(it);
        add(l, r, val);
        S.emplace_hint(it, l, r, val);
    }
    void update(const T &l, const T &r, const VAL &val) {
        update(l, r, val, [](T, T, VAL){}, [](T, T, VAL){});
    }
    template<class ADDFUNC, class DELFUNC> void insert(T l, T r, const ADDFUNC &add, const DELFUNC &del) {
        update(l, r, VAL(), add, del);
    }
    void insert(const T &l, const T &r) {
        update(l, r, VAL(), [](T, T, VAL){}, [](T, T, VAL){});
    }

    // erase [l, r)
    template<class ADDFUNC, class DELFUNC> void erase(T l, T r, const ADDFUNC &add, const DELFUNC &del) {
        auto it = S.lower_bound(Node(l, 0, VAL()));
        while (it != S.end() && it->l <= r) {
            if (it->l == r) break;
            if (it->r <= r) {
                del(it->l, it->r, it->val);
                it = S.erase(it);
            } else {
                del(it->l, r, it->val);
                Node node = *it;
                it = S.erase(it);
                it = S.emplace_hint(it, r, node.r, node.val);
            }
        }
        if (it != S.begin()) {
            it = prev(it);
            if (l < it->r) {
                if (r < it->r) {
                    it = S.emplace_hint(next(it), r, it->r, it->val);
                    it = prev(it);
                }
                del(l, min(r, it->r), it->val);
                Node node = *it;
                it = S.erase(it);
                it = S.emplace_hint(it, node.l, l, node.val);
            }
        }
    }
    void erase(const T &l, const T &r) {
        erase(l, r, [](T, T, VAL){}, [](T, T, VAL){});
    }

    // debug
    friend ostream& operator << (ostream &s, const IntervalSet &ins) {
        for (auto e : ins.S) {
            s << "([" << e.l << ", " << e.r << "): " << e.val << ") ";
        }
        return s;
    }
};


//------------------------------//
// Examples
//------------------------------//

// 第六回 アルゴリズム実技検定 M - 等しい数
void PAST_6_M() {
    long long N, Q;
    cin >> N;
    vector<long long> A(N);
    map<long long, long long> cnt;
    for (int i = 0; i < N; i++) {
        cin >> A[i];
        cnt[A[i]]++;
    }
    long long res = 0;
    for (auto [val, num] : cnt) res += num * (num - 1) / 2;

    auto add = [&](int l, int r, long long val) -> void {
        long long before = cnt[val] * (cnt[val] - 1) / 2;
        cnt[val] += r - l;
        long long after = cnt[val] * (cnt[val] - 1) / 2;
        res += after - before;
    };
    auto del = [&](int l, int r, long long val) -> void {
        long long before = cnt[val] * (cnt[val] - 1) / 2;
        cnt[val] -= r - l;
        long long after = cnt[val] * (cnt[val] - 1) / 2;
        res += after - before;
    };
    IntervalSet<int, long long> ins(A);

    cin >> Q;
    while (Q--) {
        int l, r;
        long long val;
        cin >> l >> r >> val;
        l--;
        ins.update(l, r, val, add, del);
        cout << res << '\n';
    }
}

// RUPC 2018 G - Elevator
void RUPC_2018_G() {
    using fll = array<long long, 4>;
    int N, M, Q, t, l, r;
    cin >> N >> M >> Q;
    vector<fll> qs(M + Q);
    for (int i = 0; i < M; i++) {
        cin >> t >> l >> r, l--, r--;
        qs[i] = fll({-1, t*2+1, l*2, r*2+1});
    }
    for (int i = 0; i < Q; i++) {
        cin >> t >> l >> r, l--, r--;
        qs[i+M] = fll({i, t*2, l*2, r*2});
    }
    sort(qs.begin(), qs.end(), [&](const fll &p, const fll &q) { return p[1] < q[1]; });

    IntervalSet<long long> ins;
    vector<bool> res(Q, false);
    for (auto q : qs) {
        if (q[0] == -1) {
            ins.insert(q[2], q[3]);
        } else {
            if (q[2] >= q[3] || ins.same(q[2], q[3])) res[q[0]] = true;
            else res[q[0]] = false;
        }
    }
    for (int q = 0; q < Q; ++q) cout << (res[q] ? "Yes" : "No") << endl;
}

// AtCoder ABC 255 Ex - Range Harvest Query
// modint
template<int MOD> struct Fp {
    // inner value
    long long val;
    
    // constructor
    constexpr Fp() : val(0) { }
    constexpr Fp(long long v) : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr long long get() const { return val; }
    constexpr int get_mod() const { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator + () const { return Fp(*this); }
    constexpr Fp operator - () const { return Fp(0) - Fp(*this); }
    constexpr Fp operator + (const Fp &r) const { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp &r) {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) {
        val = val * r.val % MOD;
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        val = val * u % MOD;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp pow(long long n) const {
        Fp res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const {
        Fp res(1), div(*this);
        return res / div;
    }

    // other operators
    constexpr bool operator == (const Fp &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const {
        return this->val != r.val;
    }
    constexpr Fp& operator ++ () {
        ++val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -- () {
        if (val == 0) val += MOD;
        --val;
        return *this;
    }
    constexpr Fp operator ++ (int) const {
        Fp res = *this;
        ++*this;
        return res;
    }
    constexpr Fp operator -- (int) const {
        Fp res = *this;
        --*this;
        return res;
    }
    friend constexpr istream& operator >> (istream &is, Fp<MOD> &x) {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD> &x) {
        return os << x.val;
    }
    friend constexpr Fp<MOD> pow(const Fp<MOD> &r, long long n) {
        return r.pow(n);
    }
    friend constexpr Fp<MOD> inv(const Fp<MOD> &r) {
        return r.inv();
    }
};
void ABC_255_Ex() {
    const int MOD = 998244353;
    const long long INF = 1LL << 60;
    using mint = Fp<MOD>;

    long long N, Q, D, L, R;
    cin >> N >> Q;
    mint res = 0;
    auto add = [&](long long l, long long r, long long val) -> void {};
    auto del = [&](long long l, long long r, long long val) -> void {
        res += mint(D - val) * (l + r - 1) * (r - l) / 2;
    };
    IntervalSet<long long, long long> ins;
    ins.update(0, INF, 0);
    while (Q--) {
        cin >> D >> L >> R;
        res = 0;
        ins.update(L, R+1, D, add, del);
        cout << res << '\n';
    }
}

// yukicoder No.674 n連勤
void yukicoder_674() {
    long long D, Q, A, B;
    cin >> D >> Q;
    IntervalSet<long long> ins;
    long long res = 0;
    while (Q--) {
        cin >> A >> B;
        ins.insert(A, B+1);
        auto it = ins.get(A);
        long long l = it->l, r = it->r;
        res = max(res, r - l);
        cout << res << '\n';
    }
}

// AtCoder ABC 330 E - Mex and Update
void ABC_330_E() {
    int N, Q, i, x;
    cin >> N >> Q;
    vector<int> A(N);
    map<int, int> mp;
    IntervalSet<int> ins;
    for (int i = 0; i < N; i++) {
        cin >> A[i];
        mp[A[i]]++;
        ins.insert(A[i], A[i]+1);
    }
    while (Q--) {
        cin >> i >> x;
        i--;
        mp[A[i]]--;
        if (mp[A[i]] == 0) ins.erase(A[i], A[i]+1);
        A[i] = x;
        mp[A[i]]++;
        ins.insert(A[i], A[i]+1);
        cout << ins.get_mex() << '\n';
    }
}

// 第五回 アルゴリズム実技検定 N - 旅行会社
void PAST_5_N() {
    const int IN = -2;
    const int OUT = -1;
    using tint = array<int, 3>;  // (time, type(0, 1, 2), city)
    int N, Q;
    cin >> N >> Q;
    vector<tint> events;
    for (int i = 0; i < N-1; i++) {
        int L, R;
        cin >> L >> R;
        L--;
        events.push_back({L, IN, i});
        events.push_back({R, OUT, i});
    }
    for (int i = 0; i < Q; i++) {
        int A, B;
        cin >> A >> B;
        A--, B--;
        events.push_back({A, i, B});
    }
    sort(events.begin(), events.end());

    IntervalSet<int> ins;
    vector<int> res(Q);
    for (auto [time, id, city] : events) {
        if (id == IN) {
            ins.insert(city, city+1);
        } else if (id == OUT) {
            ins.erase(city, city+1);
        } else {
            int ma = city, mi = city;
            auto it = ins.lower_bound(city);
            if (it != ins.end()) {
                if (city >= it->l) ma = max(ma, it->r), mi = min(mi, it->l);
            }
            if (it != ins.begin()) {
                it = prev(it);
                if (city <= it->r) mi = min(mi, it->l);
            }
            res[id] = ma - mi + 1;
        }
    }
    for (auto val : res) cout << val << '\n';
}

// Code Festival 2015 予選 B D - マスと駒と色塗り (for insert, lower_bound)
void code_festival_2015_B_D() {
    long long N, S, C;
    cin >> N;

    const long long INF = 1LL<<60;
    IntervalSet<long long> ins;
    ins.insert(INF, INF*2);
    while (N--) {
        cin >> S >> C;
        S--;
        long long cur = S;
        while (C > 0) {
            auto it = ins.lower_bound(cur);
            if (cur >= it->l) cur = it->r;
            else {
                long long diff = it->l - cur;
                cur += min(diff, C);
                C -= min(diff, C);
            }
        }
        ins.insert(S, cur);
        cout << cur << '\n';
    }
}

// CPSCO 2019 Session 1 E - Exclusive OR Queries
void cpsco_2019_session_1_E() {
    long long N, Q, A, L, R, X;
    cin >> N >> Q;
    set<long long> se;
    for (int i = 0; i < N; i++) {
        cin >> A;
        if (se.count(A)) se.erase(A);
        else se.insert(A);
    }
    long long res = 0, num = 0;
    auto add = [&](long long l, long long r, long long val) -> void {};
    auto del = [&](long long l, long long r, long long val) -> void {
        for (long long k = l; k < r; k++) res ^= k, num++;
    };
    IntervalSet<long long, long long> ins;
    for (auto val : se) ins.update(val, val+1, 1);
    while (Q--) {
        cin >> L >> R >> X;
        R++;
        res = 0, num = 0;
        ins.erase(L, R, add, del);
        cout << res << '\n';
        int exist = ins.covered(X);
        if ((num + exist) % 2 == 1) ins.update(X, X+1, 1);
        else ins.erase(X, X+1);
    }
}

// yukicoder No.3017 交互浴
void yukicoder_3017() {
    long long N, H;
    cin >> N;

    long long res = 0;
    auto add = [&](long long l, long long r, long long val) -> void {
        res += (r - l) * val;
    };
    auto del = [&](long long l, long long r, long long val) -> void {
        res -= (r - l) * val;
    };
    IntervalSet<long long, long long> ins;
    for (int iter = 0; iter < N; iter++) {
        cin >> H;
        if (iter % 2 == 0) ins.update(0, H, 1, add, del);
        else ins.update(0, H, 0, add, del);
        cout << res << '\n';
    }
}


int main() {
    //PAST_6_M();
    //RUPC_2018_G();
    //ABC_255_Ex();
    //yukicoder_674();
    //ABC_330_E();
    //PAST_5_N();
    //code_festival_2015_B_D();
    //cpsco_2019_session_1_E();
    yukicoder_3017();
}