//
// 抽象化した全方位木 DP
//
// verified:
//   EDPC V - Subtree
//     https://atcoder.jp/contests/dp/tasks/dp_v
//
//   TDPC N - 木
//     https://atcoder.jp/contests/tdpc/tasks/tdpc_tree
//
//   AtCoder ABC 348 E - Minimize Sum of Distances
//     https://atcoder.jp/contests/abc348/tasks/abc348_e
//

/*
    通常の木 DP において、頂点 v を根とする部分根付き木に関する再帰関数 rec(v) について、
 　　　1. res = IDENTITY
 　　　2. 頂点 v の各子頂点 v2 (その辺を e とする) に対して：res = MERGE(res, rec(v2))
 　　　3. return ADDNODE(v, res)
 　　というような更新を行うものとする。
 　　このような木 DP を全方位木 DP へと拡張する。
 */


#include <bits/stdc++.h>
using namespace std;


// re-rooting
template<class Monoid> struct ReRooting {
    using Graph = vector<vector<int>>;
    using MergeFunc = function<Monoid(Monoid, Monoid)>;
    using AddNodeFunc = function<Monoid(int, Monoid)>;
    
    // core member
    Graph G;
    Monoid IDENTITY;
    MergeFunc MERGE;
    AddNodeFunc ADDNODE;
    
    // inner data
    vector<vector<Monoid>> dp;
    vector<unordered_map<int,int>> ids;
    
    // constructor
    ReRooting() {}
    ReRooting(const Graph &g, const Monoid &identity,
              const MergeFunc &merge, const AddNodeFunc &addnode) {
        G = g;
        IDENTITY = identity;
        MERGE = merge;
        ADDNODE = addnode;
        build();
    }
    
    // re-looting dp
    Monoid rec(int v, int p) {
        Monoid res = IDENTITY;
        dp[v].assign(G[v].size(), IDENTITY);
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = G[v][i];
            ids[v][v2] = i;
            if (v2 == p) continue;
            dp[v][i] = rec(v2, v);
            res = MERGE(res, dp[v][i]);
        }
        return ADDNODE(v, res);
    }
    void rerec(int v, int p, Monoid pval) {
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = G[v][i];
            if (v2 == p) {
                dp[v][i] = pval;
                continue;
            }
        }
        vector<Monoid> left(G[v].size() + 1, IDENTITY);
        vector<Monoid> right(G[v].size() + 1, IDENTITY);
        for (int i = 0; i < G[v].size(); ++i) {
            left[i + 1] = MERGE(left[i], dp[v][i]);
            right[i + 1] = MERGE(right[i], dp[v][(int)G[v].size() - i - 1]);
        }
        for (int i = 0; i < G[v].size(); ++i) {
            int v2 = G[v][i];
            if (v2 == p) continue;
            Monoid pval2 = MERGE(left[i], right[(int)G[v].size() - i - 1]);
            rerec(v2, v, ADDNODE(v, pval2));
        }
    }
    void build() {
        dp.assign(G.size(), vector<Monoid>());
        ids.assign(G.size(), unordered_map<int,int>());
        int root = 0, nullparent = -1;
        rec(root, nullparent);
        rerec(root, nullparent, IDENTITY);
    }
    
    // getter
    Monoid get(int v) {
        Monoid res = IDENTITY;
        for (int i = 0; i < G[v].size(); ++i) {
            res = MERGE(res, dp[v][i]);
        }
        return ADDNODE(v, res);
    }
    Monoid get(int v, int w) {
        return dp[v][ids[v][w]];
    }
    
    // dump
    friend constexpr ostream& operator << (ostream &os, const ReRooting<Monoid> &rr) {
        for (int v = 0; v < rr.G.size(); ++v) {
            for (int i = 0; i < rr.G[v].size(); ++i) {
                os << v << " -> " << rr.G[v][i] << ": " << rr.dp[v][i] << endl;
            }
        }
        return os;
    }
};



//------------------------------//
// Examples
//------------------------------//

// EDPC V - Subtree
void EDPC_V() {
    int N, M;
    cin >> N >> M;
    
    using Graph = vector<vector<int>>;
    Graph G(N);
    for (int i = 0; i < N - 1; ++i) {
        int x, y;
        cin >> x >> y;
        --x, --y;
        G[x].push_back(y);
        G[y].push_back(x);
    }
    
    using Monoid = long long;
    Monoid identity = 1;
    auto merge = [&](Monoid a, Monoid b) -> Monoid { return a * b % M; };
    auto addnode = [&](int v, Monoid a) -> Monoid { return (a + 1) % M; };
    ReRooting<Monoid> rr(G, identity, merge, addnode);
    
    //cout << rr << endl;
    
    for (int v = 0; v < N; ++v) {
        cout << (rr.get(v) + M - 1) % M << endl;
    }
}

// TDPC N - 木
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

// Binomial coefficient
template<class mint> struct BiCoef {
    vector<mint> fact_, inv_, finv_;
    constexpr BiCoef() {}
    constexpr BiCoef(int n) : fact_(n, 1), inv_(n, 1), finv_(n, 1) {
        init(n);
    }
    constexpr void init(int n) {
        fact_.assign(n, 1), inv_.assign(n, 1), finv_.assign(n, 1);
        int MOD = fact_[0].get_mod();
        for(int i = 2; i < n; i++){
            fact_[i] = fact_[i-1] * i;
            inv_[i] = -inv_[MOD%i] * (MOD/i);
            finv_[i] = finv_[i-1] * inv_[i];
        }
    }
    constexpr mint com(int n, int k) const {
        if (n < k || n < 0 || k < 0) return 0;
        return fact_[n] * finv_[k] * finv_[n-k];
    }
    constexpr mint fact(int n) const {
        if (n < 0) return 0;
        return fact_[n];
    }
    constexpr mint inv(int n) const {
        if (n < 0) return 0;
        return inv_[n];
    }
    constexpr mint finv(int n) const {
        if (n < 0) return 0;
        return finv_[n];
    }
};

void TDPC_N() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;
    int N, a, b;
    cin >> N;
    vector G(N, vector<int>());
    for (int i = 0; i < N-1; i++) {
        cin >> a >> b, a--, b--;
        G[a].push_back(b), G[b].push_back(a);
    }
    BiCoef<mint> bc(N + 10);
    using Monoid = pair<int, mint>;
    Monoid identity = make_pair(0, mint(1));
    auto merge = [&](Monoid a, Monoid b) -> Monoid {
        return make_pair(a.first + b.first, a.second * b.second * bc.com(a.first + b.first, a.first));
    };
    auto addnode = [&](int v, Monoid a) -> Monoid {
        return make_pair(a.first + 1, a.second);
    };
    ReRooting<Monoid> rr(G, identity, merge, addnode);

    mint res = 0;
    for (int v = 0; v < N; v++) {
        for (auto w : G[v]) {
            if (v < w) continue;
            auto [vs, vn] = rr.get(v, w);
            auto [ws, wn] = rr.get(w, v);
            res += vn * wn * bc.com(vs + ws - 2, vs - 1);
        }
    }
    cout << res << endl;
}

// ABC 348 E - Minimize Sum of Distances
void ABC_348_E() {
    using Graph = vector<vector<int>>;
    int N;
    cin >> N;
    Graph G(N);
    vector<long long> C(N);
    for (int i = 0; i < N - 1; ++i) {
        int a, b;
        cin >> a >> b;
        --a, --b;
        G[a].push_back(b);
        G[b].push_back(a);
    }
    for (int i = 0; i < N; ++i) cin >> C[i];
    
    using Monoid = pair<long long, long long>;  // (siz, sum)
    Monoid IDENTITY = Monoid(-1, -1);
    auto MERGE = [&](Monoid a, Monoid b) {
        if (a.first == -1) return b;
        else if (b.first == -1) return a;
        else return Monoid(a.first + b.first, a.second + b.second);
    };
    auto ADDNODE = [&](int v, Monoid a) {
        Monoid res(C[v], C[v]);
        if (a.first == -1) return res;
        res.first += a.first;
        res.second += a.first + a.second;
        return res;
    };
    
    ReRooting<Monoid> rr(G, IDENTITY, MERGE, ADDNODE);
    long long res = 1LL << 62;
    for (int v = 0; v < N; ++v) {
        auto tmp = rr.get(v);
        res = min(res, tmp.second - tmp.first);
    }
    cout << res << endl;
}


int main() {
    //EDPC_V();
    TDPC_N();
    //ABC_348_E();
}