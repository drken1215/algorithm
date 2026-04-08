// code template is in https://github.com/drken1215/algorithm/blob/master/template_atcoder.cpp
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Utility
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
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl

// input
template<class T> istream& operator >> (istream &is, vector<T> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, deque<T> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, vector<vector<T>> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }

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

// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
template<class T = long long> struct Graph {
    int V, E;
    bool record_reversed_edges = false, record_edge_index = false;
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    vector<unordered_map<int, int>> id;  // id[v][w] := the index of node w in G[v]

    // constructors
    Graph(int n = 0, int m = 0, bool rre = false, bool rei = false) { 
        init(n, m, rre, rei);
    }
    void init(int n = 0, int m = 0, bool rre = false, bool rei = false) {
        V = n, E = m, record_reversed_edges = rre, record_edge_index = rei;
        list.assign(n, vector<Edge<T>>());
        if (record_reversed_edges) reversed_list.assign(n, vector<Edge<T>>());
        if (record_edge_index) id.assign(n, unordered_map<int, int>());
    }
    Graph(const Graph&) = default;
    Graph& operator = (const Graph&) = default;

    // getters
    vector<Edge<T>> &operator [] (int i) { return list[i]; }
    const vector<Edge<T>> &operator [] (int i) const { return list[i]; }
    constexpr size_t size() const { return list.size(); }
    constexpr void clear() { V = 0; list.clear(); }
    constexpr void resize(int n) { V = n; list.resize(n); }
    const vector<Edge<T>> &get_rev_edges(int i) const { 
        assert(record_reversed_edges);
        return reversed_list[i];
    }
    Edge<T> &get_edge(int u, int v) {
        assert(record_edge_index);
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u][v] >= 0 && id[u][v] < list[u].size());
        return list[u][id[u][v]];
    }
    const Edge<T> &get_edge(int u, int v) const {
        assert(record_edge_index);
        assert(u >= 0 && u < list.size() && v >= 0 && v < list.size());
        assert(id[u].count(v) && id[u].at(v) >= 0 && id[u].at(v) < list[u].size());
        return list[u][id[u].at(v)];
    }

    // add edge
    void add_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size(); 
        list[from].push_back(Edge(from, to, val));
        if (record_reversed_edges) reversed_list[to].push_back(Edge(to, from, val));
    }
    void add_bidirected_edge(int from, int to, T val = 1) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        if (record_edge_index) id[from][to] = (int)list[from].size(); 
        list[from].push_back(Edge(from, to, val));
        if (record_reversed_edges) reversed_list[from].push_back(Edge(from, to, val));
        if (from != to) {
            if (record_edge_index) id[to][from] = (int)list[to].size(); 
            list[to].push_back(Edge(to, from, val));
            if (record_reversed_edges) reversed_list[to].push_back(Edge(to, from, val));
        }
    }

    // input / output
    friend istream& operator >> (istream &is, Graph &G) {
        for (int i = 0; i < G.E; i++) {
            int u, v;
            is >> u >> v, u--, v--;
            G.add_bidirected_edge(u, v);
        }
        return is;
    }
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

// floor, ceil
template<class T> T floor(T a, T b) {
    if (a % b == 0 || a >= 0) return a / b;
    else return -((-a) / b) - 1;
}
template<class T> T ceil(T x, T y) {
    return floor(x + y - 1, y);
}

// min non-negative i such that 2^i <= n
template<class T> T floor_pow2(T n) {
    T i = 0;
    while ((T(1) << (i + 1)) <= T(n)) i++;
    return i;
}

// min non-negative i such that n <= 2^i
template<class T> T ceil_pow2(T n) {
    T i = 0;
    while ((T(1) << i) < T(n)) i++;
    return i;
}

// num of i such that (x & (1 << i)) != 0
int popcnt(int x) { return __builtin_popcount(x); }
int popcnt(unsigned int x) { return __builtin_popcount(x); }
int popcnt(long long x) { return __builtin_popcountll(x); }
int popcnt(unsigned long long x) { return __builtin_popcountll(x); }
int popcnt_mod2(int x) { return __builtin_parity(x); }
int popcnt_mod2(unsigned int x) { return __builtin_parity(x); }
int popcnt_mod2(long long x) { return __builtin_parityll(x); }
int popcnt_mod2(unsigned long long x) { return __builtin_parityll(x); }

// min non-negative i such that (x & (1 << i)) != 0
int bsf(int x) { return __builtin_ctz(x); }
int bsf(unsigned int x) { return __builtin_ctz(x); }
int bsf(long long x) { return __builtin_ctzll(x); }
int bsf(unsigned long long x) { return __builtin_ctzll(x); }

// max non-negative i such that (x & (1 << i)) != 0
int bsr(int x) { return 8 * (int)sizeof(int) - 1 - __builtin_clz(x); }
int bsr(unsigned int x) { return 8 * (int)sizeof(unsigned int) - 1 - __builtin_clz(x); }
int bsr(long long x) { return 8 * (int)sizeof(long long) - 1 - __builtin_clzll(x); }
int bsr(unsigned long long x) { return 8 * (int)sizeof(unsigned long long) - 1 - __builtin_clzll(x); }

// compression
template<class T> vector<T> compress(const vector<T> &v) {
    auto alt = v;
    sort(alt.begin(), alt.end());
    alt.erase(unique(alt.begin(), alt.end()), alt.end());
    vector<T> res(v.size());
    for (int i = 0; i < res.size(); i++) {
        res[i] = lower_bound(alt.begin(), alt.end(), v[i]) - alt.begin();
    }
    return res;
}

// Associative Array
template<class Key, class Val, uint32_t N = 20> struct FastMap {
    static constexpr uint32_t SIZE  = 1u << N;
    static constexpr uint32_t MASK  = SIZE - 1;
    static constexpr uint32_t SHIFT = 64 - N;
    static constexpr uint64_t MAGIC = 11995408973635179863ULL;
    static constexpr uint32_t get_hash(const Key& k) { return k * MAGIC >> SHIFT; }

    // inner values
    array<Key, SIZE> key;
    array<Val, SIZE> val;
    bitset<SIZE> used;

    // constructors
    FastMap() { clear(); }
    void clear() { used.reset(); }
    FastMap(const FastMap&) = default;
    FastMap& operator = (const FastMap&) = default;

    // getters
    Val &operator [] (const Key &k) {
        auto hash = get_hash(k);
        while (true) {
            if (!used[hash]) {
                used[hash] = 1;
                key[hash] = k;
                return val[hash] = Val();
            }
            if (key[hash] == k) return val[hash];
            ++hash &= MASK;
        }
    }
    const Val &operator [] (const Key &k) const {
        auto hash = get_hash(k);
        while (true) {
            if (!used[hash]) {
                used[hash] = 1;
                key[hash] = k;
                return val[hash] = Val();
            }
            if (key[hash] == k) return val[hash];
            ++hash &= MASK;
        }
    }
    Val get(const Key &k) const {
        auto hash = get_hash(k);
        while (true) {
            if (!used[hash]) return Val();
            if (key[hash] == k) return val[hash];
            ++hash &= MASK;
        }
    }
    bool count(const Key &k) const {
        auto hash = get_hash(k);
        while (true) {
            if (!used[hash]) return false;
            if (key[hash] == k) return true;
            ++hash &= MASK;
        }
    }
};

// kth root
// N < 2^64, K <= 64
uint64_t kth_root(uint64_t N, uint64_t K = 2) {
    assert(K >= 1);
    if (N <= 1 || K == 1) return N;
    if (K >= 64) return 1;
    if (N == uint64_t(-1)) --N;
    
    auto mul = [&](uint64_t x, uint64_t y) -> uint64_t {
        if (x < UINT_MAX && y < UINT_MAX) return x * y;
        if (x == uint64_t(-1) || y == uint64_t(-1)) return uint64_t(-1);
        return (x <= uint64_t(-1) / y ? x * y : uint64_t(-1));
    };
    auto power = [&](uint64_t x, uint64_t k) -> uint64_t {
        if (k == 0) return 1ULL;
        uint64_t res = 1ULL;
        while (k) {
            if (k & 1) res = mul(res, x);
            x = mul(x, x);
            k >>= 1;
        }
        return res;
    };
    
    uint64_t res;
    if (K == 2) res = sqrtl(N) - 1;
    else if (K == 3) res = cbrt(N) - 1;
    else res = pow(N, nextafter(1 / double(K), 0));
    while (power(res + 1, K) <= N) ++res;
    return res;
}

// xor128による乱数生成、周期は2^128-1
unsigned int rand_int() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty; ty = tz; tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}
int rand_int(int minv, int maxv) {
    return rand_int() % (maxv - minv + 1) + minv;
}
long long rand_ll(long long minv, long long maxv) {
    long long a = rand_int(), b = rand_int();
    return (a * (1LL<<29) + b) % (maxv - minv + 1) + minv;
}
template<class T> void shuffle(vector<T>& vec) {
    int n = vec.size();
    for (int i = n - 1; i > 0; --i) {
        int k = rand_int() % (i + 1);
        swap(vec[i], vec[k]);
    }
}

// int 128
i128 to_integer(const string &s) {
    i128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    if (s[0] == '-') res *= -1;
    return res;
}
istream& operator >> (istream &is, i128 &x) {
    string s;
    is >> s;
    x = to_integer(s);
    return is;
}
ostream& operator << (ostream &os, const i128 &x) {
    i128 ax = (x >= 0 ? x : -x);
    char buffer[128];
    char *d = end(buffer);
    do {
         --d;
        *d = "0123456789"[ax % 10];
        ax /= 10;
    } while (ax != 0);
    if (x < 0) {
        --d;
        *d = '-';
    }
    int len = end(buffer) - d;
    if (os.rdbuf()->sputn(d, len) != len) {
        os.setstate(ios_base::badbit);
    }
    return os;
}
i128 gcd(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    if (b == 0) return a;
    else return gcd(b, a % b);
}


//------------------------------//
// modint
//------------------------------//

// modint
template<int MOD = 998244353, bool PRIME = true> struct Fp {
    // inner value
    unsigned int val;
    
    // constructor
    constexpr Fp() : val(0) { }
    template<std::signed_integral T> constexpr Fp(T v) {
        long long tmp = (long long)(v % (long long)(get_umod()));
        if (tmp < 0) tmp += get_umod();
        val = (unsigned int)(tmp);
    }
    template<std::unsigned_integral T> constexpr Fp(T v) {
        val = (unsigned int)(v % get_umod());
    }
    constexpr long long get() const { return val; }
    constexpr static int get_mod() { return MOD; }
    constexpr static unsigned int get_umod() { return MOD; }
    
    // arithmetic operators
    constexpr Fp operator + () const { return Fp(*this); }
    constexpr Fp operator - () const { return Fp() - Fp(*this); }
    constexpr Fp operator + (const Fp &r) const { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp &r) const { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp &r) const { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp &r) const { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp &r) {
        val += r.val;
        if (val >= get_umod()) val -= get_umod();
        return *this;
    }
    constexpr Fp& operator -= (const Fp &r) {
        val -= r.val;
        if (val >= get_umod()) val += get_umod();
        return *this;
    }
    constexpr Fp& operator *= (const Fp &r) {
        unsigned long long tmp = val;
        tmp *= r.val;
        val = (unsigned int)(tmp % get_umod());
        return *this;
    }
    constexpr Fp& operator /= (const Fp &r) {
        return *this = *this * r.inv(); 
    }
    constexpr Fp pow(long long n) const {
        assert(n >= 0);
        Fp res(1), mul(*this);
        while (n) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    constexpr Fp inv() const {
        assert(val);
        if (PRIME) {
            return pow(get_umod() - 2);
        } else {
            assert(gcd(val, get_umod()) == 1);
            long long m = get_umod(), a = val, b = m, u = 1, v = 0;
            while (b > 0) {
                auto t = a / b;
                a -= t * b, swap(a, b);
                u -= t * v, swap(u, v);
            }
            return Fp(u);
        }
    }

    // other operators
    constexpr bool operator == (const Fp &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp &r) const {
        return this->val != r.val;
    }
    constexpr bool operator < (const Fp &r) const {
        return this->val < r.val;
    }
    constexpr bool operator > (const Fp &r) const {
        return this->val > r.val;
    }
    constexpr bool operator <= (const Fp &r) const {
        return this->val <= r.val;
    }
    constexpr bool operator >= (const Fp &r) const {
        return this->val >= r.val;
    }
    constexpr Fp& operator ++ () {
        ++val;
        if (val == get_umod()) val = 0;
        return *this;
    }
    constexpr Fp& operator -- () {
        if (val == 0) val = get_umod();
        --val;
        return *this;
    }
    constexpr Fp operator ++ (int) {
        Fp res = *this;
        ++*this;
        return res;
    }
    constexpr Fp operator -- (int) {
        Fp res = *this;
        --*this;
        return res;
    }
    friend constexpr istream& operator >> (istream &is, Fp &x) {
        long long tmp = 1;
        is >> tmp;
        tmp = tmp % (long long)(get_umod());
        if (tmp < 0) tmp += get_umod();
        x.val = (unsigned int)(tmp);
        return is;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp &x) {
        return os << x.val;
    }
    friend constexpr Fp pow(const Fp &r, long long n) {
        return r.pow(n);
    }
    friend constexpr Fp inv(const Fp &r) {
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

// dynamic modint
struct DynamicModint {
    using mint = DynamicModint;
    
    // static menber
    static int MOD;
    
    // inner value
    unsigned int val;
    
    // constructor
    DynamicModint() : val(0) { }
    template<std::signed_integral T> DynamicModint(T v) {
        long long tmp = (long long)(v % (long long)(get_umod()));
        if (tmp < 0) tmp += get_umod();
        val = (unsigned int)(tmp);
    }
    template<std::unsigned_integral T> DynamicModint(T v) {
        val = (unsigned int)(v % get_umod());
    }
    long long get() const { return val; }
    static int get_mod() { return MOD; }
    static unsigned int get_umod() { return MOD; }
    static void set_mod(int mod) { MOD = mod; }
    
    // arithmetic operators
    mint operator + () const { return mint(*this); }
    mint operator - () const { return mint() - mint(*this); }
    mint operator + (const mint &r) const { return mint(*this) += r; }
    mint operator - (const mint &r) const { return mint(*this) -= r; }
    mint operator * (const mint &r) const { return mint(*this) *= r; }
    mint operator / (const mint &r) const { return mint(*this) /= r; }
    mint& operator += (const mint &r) {
        val += r.val;
        if (val >= get_umod()) val -= get_umod();
        return *this;
    }
    mint& operator -= (const mint &r) {
        val -= r.val;
        if (val >= get_umod()) val += get_umod();
        return *this;
    }
    mint& operator *= (const mint &r) {
        unsigned long long tmp = val;
        tmp *= r.val;
        val = (unsigned int)(tmp % get_umod());
        return *this;
    }
    mint& operator /= (const mint &r) {
        return *this = *this * r.inv(); 
    }
    mint pow(long long n) const {
        assert(n >= 0);
        mint res(1), mul(*this);
        while (n) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    mint inv() const {
        assert(val);
        assert(gcd(val, get_umod()) == 1);
        long long m = get_umod(), a = val, b = m, u = 1, v = 0;
        while (b > 0) {
            auto t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        return mint(u);
    }

    // other operators
    bool operator == (const mint &r) const {
        return this->val == r.val;
    }
    bool operator != (const mint &r) const {
        return this->val != r.val;
    }
    bool operator < (const mint &r) const {
        return this->val < r.val;
    }
    bool operator > (const mint &r) const {
        return this->val > r.val;
    }
    bool operator <= (const mint &r) const {
        return this->val <= r.val;
    }
    bool operator >= (const mint &r) const {
        return this->val >= r.val;
    }
    mint& operator ++ () {
        ++val;
        if (val == get_umod()) val = 0;
        return *this;
    }
    mint& operator -- () {
        if (val == 0) val = get_umod();
        --val;
        return *this;
    }
    mint operator ++ (int) {
        mint res = *this;
        ++*this;
        return res;
    }
    mint operator -- (int) {
        mint res = *this;
        --*this;
        return res;
    }
    friend istream& operator >> (istream &is, mint &x) {
        long long tmp = 1;
        is >> tmp;
        tmp = tmp % (long long)(get_umod());
        if (tmp < 0) tmp += get_umod();
        x.val = (unsigned int)(tmp);
        return is;
    }
    friend ostream& operator << (ostream &os, const mint &x) {
        return os << x.val;
    }
    friend mint pow(const mint &r, long long n) {
        return r.pow(n);
    }
    friend mint inv(const mint &r) {
        return r.inv();
    }
};
int DynamicModint::MOD;


//------------------------------//
// Number Theory
//------------------------------//

// mod pow
template<class T_VAL, class T_MOD>
constexpr T_VAL mod_pow(T_VAL a, T_VAL n, T_MOD m) {
    T_VAL res = 1;
    while (n > 0) {
        if (n % 2 == 1) res = res * a % m;
        a = a * a % m;
        n >>= 1;
    }
    return res;
}

// isprime[n] := is n prime?
// mebius[n] := mebius value of n
// min_factor[n] := the min prime-factor of n
// euler[n] := euler function value of n
struct Eratos {
    int N;
    vector<int> primes;
    vector<bool> isprime;
    vector<int> mebius, min_factor, euler;

    // constructor, getter
    Eratos(int n) : N(n), 
                    primes(),
                    isprime(n, true),
                    mebius(n, 1),
                    min_factor(n, -1),
                    euler(n) {
        isprime[0] = isprime[1] = false;
        min_factor[0] = 0, min_factor[1] = 1;
        for (int i = 1; i < N; i++) euler[i] = i;
        for (int i = 2; i < N; ++i) {
            if (!isprime[i]) continue;
            primes.push_back(i);
            mebius[i] = -1;
            min_factor[i] = i;
            euler[i] = i - 1;
            for (int j = i*2; j < N; j += i) {
                isprime[j] = false;
                if ((j / i) % i == 0) mebius[j] = 0;
                else mebius[j] = -mebius[j];
                if (min_factor[j] == -1) min_factor[j] = i;
                euler[j] /= i, euler[j] *= i - 1;
            }
        }
    }

    // prime factorization
    vector<pair<int,int>> prime_factors(int n) {
        assert(n >= 1 && n < N);
        vector<pair<int,int>> res;
        while (n != 1) {
            int prime = min_factor[n];
            int exp = 0;
            while (min_factor[n] == prime) {
                ++exp;
                n /= prime;
            }
            res.push_back(make_pair(prime, exp));
        }
        return res;
    }

    // enumerate divisors
    vector<int> divisors(int n) {
        assert(n >= 1 && n < N);
        vector<int> res({1});
        auto pf = prime_factors(n);
        for (auto p : pf) {
            int n = (int)res.size();
            for (int i = 0; i < n; ++i) {
                int v = 1;
                for (int j = 0; j < p.second; ++j) {
                    v *= p.first;
                    res.push_back(res[i] * v);
                }
            }
        }
        return res;
    }
};

// montgomery modint (MOD < 2^62, MOD: odd prime number)
struct MontgomeryModInt64 {
    using mint = MontgomeryModInt64;
    using u64 = uint64_t;
    using u128 = __uint128_t;
    
    // static menber
    static u64 MOD;
    static u64 INV_MOD;  // INV_MOD * MOD ≡ 1 (mod 2^64)
    static u64 T128;  // 2^128 (mod MOD)
    
    // inner value
    u64 val;
    
    // constructor
    MontgomeryModInt64() : val(0) { }
    MontgomeryModInt64(long long v) : val(reduce((u128(v) + MOD) * T128)) { }
    u64 get() const {
        u64 res = reduce(val);
        return res >= MOD ? res - MOD : res;
    }
    
    // mod getter and setter
    static u64 get_mod() { return MOD; }
    static void set_mod(u64 mod) {
        assert(mod < (1LL << 62));
        assert((mod & 1));
        MOD = mod;
        T128 = -u128(mod) % mod;
        INV_MOD = get_inv_mod();
    }
    static u64 get_inv_mod() {
        u64 res = MOD;
        for (int i = 0; i < 5; ++i) res *= 2 - MOD * res;
        return res;
    }
    static u64 reduce(const u128 &v) {
        return (v + u128(u64(v) * u64(-INV_MOD)) * MOD) >> 64;
    }
    
    // arithmetic operators
    mint operator + () const { return mint(*this); }
    mint operator - () const { return mint() - mint(*this); }
    mint operator + (const mint &r) const { return mint(*this) += r; }
    mint operator - (const mint &r) const { return mint(*this) -= r; }
    mint operator * (const mint &r) const { return mint(*this) *= r; }
    mint operator / (const mint &r) const { return mint(*this) /= r; }
    mint& operator += (const mint &r) {
        if ((val += r.val) >= 2 * MOD) val -= 2 * MOD;
        return *this;
    }
    mint& operator -= (const mint &r) {
        if ((val += 2 * MOD - r.val) >= 2 * MOD) val -= 2 * MOD;
        return *this;
    }
    mint& operator *= (const mint &r) {
        val = reduce(u128(val) * r.val);
        return *this;
    }
    mint& operator /= (const mint &r) {
        *this *= r.inv();
        return *this;
    }
    mint pow(u128 n) const {
        mint res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    mint inv() const { 
        return pow(MOD - 2); 
    }

    // other operators
    bool operator == (const mint &r) const {
        return (val >= MOD ? val - MOD : val) == (r.val >= MOD ? r.val - MOD : r.val);
    }
    bool operator != (const mint &r) const {
        return (val >= MOD ? val - MOD : val) != (r.val >= MOD ? r.val - MOD : r.val);
    }
    mint& operator ++ () {
        ++val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    mint& operator -- () {
        if (val == 0) val += MOD;
        --val;
        return *this;
    }
    mint operator ++ (int) {
        mint res = *this;
        ++*this;
        return res;
    }
    mint operator -- (int) {
        mint res = *this;
        --*this;
        return res;
    }
    friend istream& operator >> (istream &is, mint &x) {
        long long t;
        is >> t;
        x = mint(t);
        return is;
    }
    friend ostream& operator << (ostream &os, const mint &x) {
        return os << x.get();
    }
    friend mint pow(const mint &r, long long n) {
        return r.pow(n);
    }
    friend mint inv(const mint &r) {
        return r.inv();
    }
};

typename MontgomeryModInt64::u64
MontgomeryModInt64::MOD, MontgomeryModInt64::INV_MOD, MontgomeryModInt64::T128;

// Miller-Rabin
bool MillerRabin(long long N, const vector<long long> &A) {
    assert(N % 2 == 1);
    assert(N < (1LL<<62));
    using mint = MontgomeryModInt64;
    mint::set_mod(N);
    
    long long s = 0, d = N - 1;
    while (d % 2 == 0) {
        ++s;
        d >>= 1;
    }
    for (auto a : A) {
        if (N <= a) return true;
        mint x = mint(a).pow(d);
        if (x != 1) {
            long long t;
            for (t = 0; t < s; ++t) {
                if (x == N - 1) break;
                x *= x;
            }
            if (t == s) return false;
        }
    }
    return true;
}

bool is_prime(long long N) {
    if (N <= 1) return false;
    else if (N == 2) return true;
    else if (N % 2 == 0) return false;
    else if (N < 4759123141LL)
        return MillerRabin(N, {2, 7, 61});
    else
        return MillerRabin(N, {2, 325, 9375, 28178, 450775, 9780504, 1795265022});
}

// Pollard's Rho
unsigned int xor_shift_rng() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty, ty = tz, tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}

long long pollard(long long N) {
    if (N % 2 == 0) return 2;
    if (is_prime(N)) return N;
    
    assert(N < (1LL<<62));
    using mint = MontgomeryModInt64;
    mint::set_mod(N);
    
    long long step = 0;
    while (true) {
        mint r = xor_shift_rng();  // random r
        auto f = [&](mint x) -> mint { return x * x + r; };
        mint x = ++step, y = f(x);
        while (true) {
            long long p = gcd((y - x).get(), N);
            if (p == 0 || p == N) break;
            if (p != 1) return p;
            x = f(x);
            y = f(f(y));
        }
    }
}

vector<long long> pollard_prime_factorize(long long N) {
    if (N == 1) return {};
    long long p = pollard(N);
    if (p == N) return {p};
    vector<long long> left = pollard_prime_factorize(p);
    vector<long long> right = pollard_prime_factorize(N / p);
    if (left.size() > right.size()) swap(left, right);
    left.insert(left.end(), right.begin(), right.end());
    sort(left.begin(), left.end());
    return left;
}

vector<pair<long long, long long>> prime_factorize(long long N) {
    vector<pair<long long, long long>> res;
    const auto &prs = pollard_prime_factorize(N);
    long long prev = -1, num = 0;
    for (const auto &pr : prs) {
        if (pr == prev) ++num;
        else {
            if (prev != -1) res.emplace_back(prev, num);
            prev = pr, num = 1;
        }
    }
    if (prev != -1) res.emplace_back(prev, num);
    return res;
}

// ax + by = gcd(a, b)
// return: gcd(a, b)
template<class T_VAL> T_VAL ext_gcd(T_VAL a, T_VAL b, T_VAL &x, T_VAL &y) {
    T_VAL xsign = 1, ysign = 1;
    if (a < 0) {
        a = -a, xsign *= -1;
    }
    if (b < 0) {
        b = -b, ysign *= -1;
    }
    
    if (b == 0) {
        x = xsign, y = 0;
        return a;
    }
    T_VAL d = ext_gcd(b, a % b, y, x);
    y -= a / b * x;
    x *= xsign, y *= ysign;
    return d;
}

// calc order
template<class T_VAL, class T_MOD>
constexpr T_VAL calc_order(T_VAL a, T_MOD m) {
    assert(gcd(a, m) == 1);
    T_VAL euler = m;
    for (const auto &[p, exp] : prime_factorize(m)) euler -= euler / p;
    T_VAL res = euler;
    for (const auto &[p, exp] : prime_factorize(euler)) {
        while (res % p == 0 && mod_pow(a, res / p, m) == 1) res /= p;
    }
    return res;
}

// Discrete Logarithm (in O(√m))
// find non-negative x s.t. a^x = b (mod m)
// return (c, d) where solution is x = c (mod d)
template<class T_VAL, class T_MOD>
constexpr pair<T_VAL, T_MOD> mod_log(T_VAL a, T_VAL b, T_MOD m) {
    if (m == 1) return {0, 1};
    if (a == 0 && b == 0) return {1, 1};
    T_VAL m1 = m, pw = 1, div = gcd(a, m1), res = 0, x, y;
    for (; (div = gcd(a, m1)) > 1; res++, m1 /= div, pw = pw * a % m) {
        if (pw == b) return {res, 0};  // aperiodic solution
    }
    auto g = ext_gcd(pw, m, x, y);
    if (b % g > 0) return {-1, 0};  // no solution
    b = x * (b / g) % m;
    if (b < 0) b += m;
    T_VAL order = calc_order(a, m1), sq = kth_root(order, 2) + 1, an = mod_pow(a, sq, m1);
    FastMap<T_VAL, T_VAL> half;
    pw = an;
    for (T_VAL p = 1; p <= sq; p++, pw = pw * an % m1) half[pw] = p;
    pw = b % m1;
    for (T_VAL q = 0; q <= sq; q++, pw = pw * a % m1) {
        if (half.count(pw)) {
            res += (sq * half[pw] - q) % order;
            return {res, order};  // periodic solution
        }
    }
    return {-1, 0};  // no solution
}


//------------------------------//
// NTT
//------------------------------//

// calc primitive root
constexpr int calc_primitive_root(long long m) {
    if (m == 1) return -1;
    if (m == 2) return 1;
    if (m == 998244353) return 3;
    if (m == 167772161) return 3;
    if (m == 469762049) return 3;
    if (m == 754974721) return 11;
    if (m == 645922817) return 3;
    if (m == 897581057) return 3;

    auto mod_pow = [&](long long a, long long n, long long m) {
        long long res = 1;
        while (n > 0) {
            if (n % 2 == 1) res = res * a % m;
            a = a * a % m;
            n >>= 1;
        }
        return res;
    };
    long long divs[20] = {};
    divs[0] = 2;
    long long cnt = 1;
    long long x = (m - 1) / 2;
    while (x % 2 == 0) x /= 2;
    for (long long i = 3; i * i <= x; i += 2) {
        if (x % i == 0) {
            divs[cnt++] = i;
            while (x % i == 0) x /= i;
        }
    }
    if (x > 1) divs[cnt++] = x;
    for (long long g = 2; ; g++) {
        bool ok = true;
        for (int i = 0; i < cnt; i++) {
            if (mod_pow(g, (m - 1) / divs[i], m) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
}

// NTT setup
template<class mint, int MOD = mint::get_mod(), int g = calc_primitive_root(mint::get_mod())>
struct ntt_setup {
    static constexpr int bsf_constexpr(unsigned int x) {
        int i = 0;
        while (!(x & (1 << i))) i++;
        return i;
    };

    static constexpr int rank = bsf_constexpr(MOD - 1);
    array<mint, rank + 1> root, iroot;  // root[i]^(2^i) = 1, root[i] * iroot[i] = 1
    array<mint, max(0, rank - 1)> rate2, irate2;
    array<mint, max(0, rank - 2)> rate3, irate3;

    ntt_setup() {
        root[rank] = mint(g).pow((MOD - 1) >> rank);
        iroot[rank] = root[rank].inv();
        for (int i = rank - 1; i >= 0; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }
        mint prod = 1, iprod = 1;
        for (int i = 0; i < rank - 1; i++) {
            rate2[i] = root[i + 2] * prod;
            irate2[i] = iroot[i + 2] * iprod;
            prod *= iroot[i + 2];
            iprod *= root[i + 2];
        }
        prod = 1, iprod = 1;
        for (int i = 0; i < rank - 2; i++) {
            rate3[i] = root[i + 3] * prod;
            irate3[i] = iroot[i + 3] * iprod;
            prod *= iroot[i + 3];
            iprod *= root[i + 3];
        }
    }
};

// NTT transformation
template<class mint, int MOD = mint::get_mod()> 
void ntt_trans(vector<mint> &v) {
    int n = (int)v.size();
    int h = 0;
    while ((1U << h) < (unsigned int)(n)) h++;
    static const ntt_setup<mint> setup;

    int len = 0;
    while (len < h) {
        if (h - len == 1) {
            int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); s++) {
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto l = v[i + offset];
                    auto r = v[i + offset + p] * rot;
                    v[i + offset] = l + r;
                    v[i + offset + p] = l - r;
                }
                if (s + 1 != (1 << len)) {
                    rot *= setup.rate2[setup.bsf_constexpr(~(unsigned int)(s))];
                }
            }
            len++;
        } else {
            int p = 1 << (h - len - 2);
            mint rot = 1, imag = setup.root[2];
            for (int s = 0; s < (1 << len); s++) {
                mint rot2 = rot * rot, rot3 = rot2 * rot;
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto mod2 = 1ULL * MOD * MOD;
                    auto a0 = 1ULL * v[i + offset].val;
                    auto a1 = 1ULL * v[i + offset + p].val * rot.val;
                    auto a2 = 1ULL * v[i + offset + p * 2].val * rot2.val;
                    auto a3 = 1ULL * v[i + offset + p * 3].val * rot3.val;
                    auto tmp = 1ULL * mint(a1 + mod2 - a3).val * imag.val;
                    auto na2 = mod2 - a2;
                    v[i + offset] = a0 + a2 + a1 + a3;
                    v[i + offset + p] = a0 + a2 + (mod2 * 2 - (a1 + a3));
                    v[i + offset + p * 2] = a0 + na2 + tmp;
                    v[i + offset + p * 3] = a0 + na2 + (mod2 - tmp);
                }
                if (s + 1 != (1 << len)) {
                    rot *= setup.rate3[setup.bsf_constexpr(~(unsigned int)(s))];
                }
            }
            len += 2;
        }
    }
}

// NTT inv-transformation
template<class mint, int MOD = mint::get_mod()> 
void ntt_trans_inv(vector<mint> &v) {
    int n = (int)v.size();
    int h = 0;
    while ((1U << h) < (unsigned int)(n)) h++;
    static const ntt_setup<mint> setup;

    int len = h;
    while (len) {
        if (len == 1) {
            int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); s++) {
                int offset = s << (h - len + 1);
                for (int i = 0; i < p; i++) {
                    auto l = v[i + offset];
                    auto r = v[i + offset + p];
                    v[i + offset] = l + r;
                    v[i + offset + p] = (unsigned long long)((long long)(MOD) + l.val - r.val) * irot.val;
                }
                if (s + 1 != (1 << (len - 1))) {
                    irot *= setup.irate2[setup.bsf_constexpr(~(unsigned int)(s))];
                }
            }
            len--;
        } else {
            int p = 1 << (h - len);
            mint irot = 1, iimag = setup.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); s++) {
                mint irot2 = irot * irot, irot3 = irot2 * irot;
                int offset = s << (h - len + 2);
                for (int i = 0; i < p; i++) {
                    auto a0 = 1ULL * v[i + offset].val;
                    auto a1 = 1ULL * v[i + offset + p].val;
                    auto a2 = 1ULL * v[i + offset + p * 2].val;
                    auto a3 = 1ULL * v[i + offset + p * 3].val;
                    auto tmp = 1ULL * mint((MOD + a2 - a3) * iimag.val).val;
                    v[i + offset] = a0 + a1 + a2 + a3;
                    v[i + offset + p] = (a0 + (MOD - a1) + tmp) * irot.val;
                    v[i + offset + p * 2] = (a0 + a1 + (MOD - a2) + (MOD - a3)) * irot2.val;
                    v[i + offset + p * 3] = (a0 + (MOD - a1) + (MOD - tmp)) * irot3.val;
                }
                if (s + 1 != (1 << (len - 2))) {
                    irot *= setup.irate3[setup.bsf_constexpr(~(unsigned int)(s))];
                }
            }
            len -= 2;
        }
    }
    mint in = mint(n).inv();
    for (int i = 0; i < n; i++) v[i] *= in;
}

// naive convolution
template<class VEC> VEC convolution_naive(const VEC &a, const VEC &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    VEC res(n + m - 1);
    if (n < m) {
        for (int j = 0; j < m; j++) for (int i = 0; i < n; i++) res[i + j] += a[i] * b[j];
    } else {
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) res[i + j] += a[i] * b[j];
    }
    return res;
}

// ntt convolution
template<class mint>
vector<mint> convolution_ntt(vector<mint> a, vector<mint> b) {
    int MOD = mint::get_mod();
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    int z = (int)bit_ceil((unsigned int)(n + m - 1));
    assert((MOD - 1) % z == 0);
    a.resize(z), b.resize(z);
    ntt_trans(a), ntt_trans(b);
    for (int i = 0; i < z; i++) a[i] *= b[i];
    ntt_trans_inv(a);
    a.resize(n + m - 1);
    return a;
}

// convolution long long (if u64 is necessary, use convolution_ull)
template<class VEC> VEC convolution_ll(const VEC &a, const VEC &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return VEC();
    if (min(n, m) <= 60) return convolution_naive(a, b);

    static constexpr int MOD0 = 754974721;  // 2^24
    static constexpr int MOD1 = 167772161;  // 2^25
    static constexpr int MOD2 = 469762049;  // 2^26
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    static const mint1 imod0 = 95869806; // modinv(MOD0, MOD1);
    static const mint2 imod1 = 104391568; // modinv(MOD1, MOD2);
    static const mint2 imod01 = 187290749; // imod1 / MOD0;

    vector<mint0> a0(n, 0), b0(m, 0);
    vector<mint1> a1(n, 0), b1(m, 0);
    vector<mint2> a2(n, 0), b2(m, 0);
    for (int i = 0; i < n; ++i) a0[i] = a[i], a1[i] = a[i], a2[i] = a[i];
    for (int i = 0; i < m; ++i) b0[i] = b[i], b1[i] = b[i], b2[i] = b[i];
    auto c0 = convolution_ntt(std::move(a0), std::move(b0));
    auto c1 = convolution_ntt(std::move(a1), std::move(b1));
    auto c2 = convolution_ntt(std::move(a2), std::move(b2));

    VEC res(n + m - 1);
    long long mod0 = MOD0, mod01 = mod0 * MOD1;
    for (int i = 0; i < n + m - 1; ++i) {
        unsigned int y0 = c0[i].val;
        unsigned int y1 = (imod0 * (c1[i] - mint1(y0))).val;
        unsigned int y2 = (imod01 * (c2[i] - mint2(y0)) - imod1 * y1).val;
        res[i] = mod01 * y2 + mod0 * y1 + y0;
    }
    return res;
}

// convolution in general mod
template<class mint>
vector<mint> convolution_general_mod(const vector<mint> &a, const vector<mint> &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    if (min(n, m) <= 60) return convolution_naive(a, b);
    if constexpr (std::is_same_v<mint, Fp<998244353>>) return convolution_ntt(a, b);

    static constexpr int MOD0 = 754974721;  // 2^24
    static constexpr int MOD1 = 167772161;  // 2^25
    static constexpr int MOD2 = 469762049;  // 2^26
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    static const mint1 imod0 = 95869806; // modinv(MOD0, MOD1);
    static const mint2 imod1 = 104391568; // modinv(MOD1, MOD2);
    static const mint2 imod01 = 187290749; // imod1 / MOD0;

    vector<mint0> a0(n, 0), b0(m, 0);
    vector<mint1> a1(n, 0), b1(m, 0);
    vector<mint2> a2(n, 0), b2(m, 0);
    for (int i = 0; i < n; ++i) a0[i] = a[i].val, a1[i] = a[i].val, a2[i] = a[i].val;
    for (int i = 0; i < m; ++i) b0[i] = b[i].val, b1[i] = b[i].val, b2[i] = b[i].val;
    auto c0 = convolution_ntt(std::move(a0), std::move(b0));
    auto c1 = convolution_ntt(std::move(a1), std::move(b1));
    auto c2 = convolution_ntt(std::move(a2), std::move(b2));

    vector<mint> res(n + m - 1);
    mint mod0 = MOD0, mod01 = mod0 * MOD1;
    for (int i = 0; i < n + m - 1; ++i) {
        unsigned int y0 = c0[i].val;
        unsigned int y1 = (imod0 * (c1[i] - mint1(y0))).val;
        unsigned int y2 = (imod01 * (c2[i] - mint2(y0)) - imod1 * y1).val;
        res[i] = mod01 * y2 + mod0 * y1 + y0;
    }
    return res;
}

// convolution overall
template<class T> vector<T> convolution(const vector<T> &a, const vector<T> &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    if (min(n, m) <= 60) return convolution_naive(a, b);
    if constexpr (std::is_same_v<T, Fp<998244353>>) return convolution_ntt(a, b);
    else if constexpr (std::is_integral_v<T>) return convolution_ll(a, b);
    else return convolution_general_mod(a, b);
}


//------------------------------//
// FPS
//------------------------------//

// mod sqrt
template<class T_VAL, class T_MOD>
T_VAL mod_sqrt(T_VAL a, T_MOD p) {
    a = safe_mod(a, p);
    if (a <= 1) return a;
    using mint = DynamicModint;
    mint::set_mod(p);
    if (mint(a).pow((p - 1) >> 1) != 1) return T_VAL(-1);
    mint b = 1, one = 1;
    while (b.pow((p - 1) >> 1) == 1) b++;
    T_VAL m = p - 1, e = 0;
    while (m % 2 == 0) m >>= 1, e++;
    mint x = mint(a).pow((m - 1) >> 1);
    mint y = mint(a) * x * x;
    x *= a;
    mint z = mint(b).pow(m);
    while (y != 1) {
        T_VAL j = 0;
        mint t = y;
        while (t != one) {
            j++;
            t *= t;
        }
        z = z.pow(T_VAL(1) << (e - j - 1));
        x *= z, z *= z, y *= z;
        e = j;
    }
    T_VAL res = x.val;
    if (res * 2 > p) res = p - res;
    return res;
}

// Formal Power Series
template<class mint> struct FPS : vector<mint> {
    static const int SPARSE_BOARDER = 60;
    using vector<mint>::vector;
 
    // constructor
    constexpr FPS(const vector<mint> &r) : vector<mint>(r) {}
 
    // core operator
    constexpr FPS pre(int siz) const {
        return FPS(begin(*this), begin(*this) + min((int)this->size(), siz));
    }
    constexpr FPS rev() const {
        FPS res = *this;
        reverse(begin(res), end(res));
        return res;
    }
    constexpr FPS& normalize() {
        while (!this->empty() && this->back() == 0) this->pop_back();
        return *this;
    }
    constexpr mint eval(const mint &v) const {
        mint res = 0;
        for (int i = (int)this->size()-1; i >= 0; --i) {
            res *= v;
            res += (*this)[i];
        }
        return res;
    }
    constexpr int count_terms() const {
        int res = 0;
        for (int i = 0; i < (int)this->size(); i++) if ((*this)[i] != mint(0)) res++;
        return res;
    }
 
    // basic operator
    constexpr FPS operator - () const noexcept {
        FPS res = (*this);
        for (int i = 0; i < (int)res.size(); ++i) res[i] = -res[i];
        return res;
    }
    constexpr FPS operator + (const mint &v) const { return FPS(*this) += v; }
    constexpr FPS operator + (const FPS &r) const { return FPS(*this) += r; }
    constexpr FPS operator - (const mint &v) const { return FPS(*this) -= v; }
    constexpr FPS operator - (const FPS &r) const { return FPS(*this) -= r; }
    constexpr FPS operator * (const mint &v) const { return FPS(*this) *= v; }
    constexpr FPS operator * (const FPS &r) const { return FPS(*this) *= r; }
    constexpr FPS operator / (const mint &v) const { return FPS(*this) /= v; }
    constexpr FPS operator / (const FPS &r) const { return FPS(*this) /= r; }
    constexpr FPS operator % (const FPS &r) const { return FPS(*this) %= r; }
    constexpr FPS operator << (int x) const { return FPS(*this) <<= x; }
    constexpr FPS operator >> (int x) const { return FPS(*this) >>= x; }
    constexpr FPS& operator += (const mint &v) {
        if (this->empty()) this->reserve(1), this->resize(1);
        (*this)[0] += v;
        return *this;
    }
    constexpr FPS& operator += (const FPS &r) {
        if (r.size() > this->size()) this->reserve(r.size()), this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return this->normalize();
    }
    constexpr FPS& operator -= (const mint &v) {
        if (this->empty()) this->reserve(1), this->resize(1);
        (*this)[0] -= v;
        return *this;
    }
    constexpr FPS& operator -= (const FPS &r) {
        if (r.size() > this->size()) this->reserve(r.size()), this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return this->normalize();
    }
    constexpr FPS& operator *= (const mint &v) {
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= v;
        return *this;
    }
    constexpr FPS& operator *= (const FPS &r) {
        return *this = convolution((*this), r);
    }
    constexpr FPS& divide_by_modint(const mint &v) {
        assert(v != 0);
        mint iv = v.inv();
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= iv;
        return *this;
    }
    constexpr FPS& divide_by_integer(const mint &v) {
        assert(v != 0);
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] /= v;
        return *this;
    }
    constexpr FPS& operator /= (const mint &v) {
        assert(v != 0);
        if constexpr (std::is_integral_v<mint>) return divide_by_integer(v);
        else return divide_by_modint(v);
    }
    
    // division, r must be normalized (r.back() must not be 0)
    constexpr FPS& operator /= (const FPS &r) {
        assert(!r.empty());
        assert(r.back() != 0);
        this->normalize();
        if (this->size() < r.size()) {
            this->clear();
            return *this;
        }
        int need = (int)this->size() - (int)r.size() + 1;
        *this = (rev().pre(need) * r.rev().inv(need)).pre(need).rev();
        return *this;
    }
    constexpr FPS& operator %= (const FPS &r) {
        assert(!r.empty());
        assert(r.back() != 0);
        this->normalize();
        FPS q = (*this) / r;
        return *this -= q * r;
    }
    constexpr FPS& operator <<= (int x) {
        FPS res(x, 0);
        res.insert(res.end(), begin(*this), end(*this));
        return *this = res;
    }
    constexpr FPS& operator >>= (int x) {
        FPS res;
        res.insert(res.end(), begin(*this) + x, end(*this));
        return *this = res;
    }

    // advanced operation
    // df/dx
    constexpr FPS diff() const {
        int n = (int)this->size();
        if (n <= 0) return FPS();
        FPS res(n-1);
        for (int i = 1; i < n; ++i) res[i-1] = (*this)[i] * i;
        return res;
    }
    
    // \int f dx
    constexpr FPS integral() const {
        int n = (int)this->size();
        FPS res(n+1, 0);
        for (int i = 0; i < n; ++i) res[i+1] = (*this)[i] / (i+1);
        return res;
    }
    
    // inv(f), f[0] must not be 0
    constexpr FPS inv(int deg = -1) const {
        if (count_terms() <= SPARSE_BOARDER) return inv_sparse(deg);
        if constexpr (std::is_same_v<mint, Fp<998244353>>) return inv_ntt_friendly(deg);
        assert(this->size() >= 1 && (*this)[0] != 0);
        if (deg < 0) deg = (int)this->size();
        FPS res({mint(1) / (*this)[0]});
        for (int d = 1; d < deg; d <<= 1) {
            res = (res + res - res * res * pre(d << 1)).pre(d << 1);
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS inv_ntt_friendly(int deg = -1) const {
        assert(this->size() >= 1 && (*this)[0] != 0);
        if (deg < 0) deg = (int)this->size();
        FPS res(deg);
        res[0] = mint(1) / (*this)[0];
        for (int d = 1; d < deg; d <<= 1) {
            FPS g(d * 2), h(d * 2);
            mint iv = mint(d * 2).inv();
            for (int i = 0; i < min((int)this->size(), d * 2); i++) g[i] = (*this)[i];
            for (int i = 0; i < d; i++) h[i] = res[i];
            ntt_trans(g), ntt_trans(h);
            for (int i = 0; i < d * 2; i++) g[i] *= h[i];
            ntt_trans_inv(g);
            for (int i = 0; i < d; i++) g[i] = 0;
            ntt_trans(g);
            for (int i = 0; i < d * 2; i++) g[i] *= h[i];
            ntt_trans_inv(g);
            for (int i = d; i < min(deg, d * 2); i++) res[i] = -g[i];
        }
        return res.pre(deg);
    }
    constexpr FPS inv_sparse(int deg = -1) const {
        assert(this->size() >= 1 && (*this)[0] != 0);
        if (deg < 0) deg = (int)this->size();
        vector<pair<int, mint>> dat;
        for (int i = 1; i < (int)this->size(); i++) if ((*this)[i] != mint(0)) {
            dat.emplace_back(i, (*this)[i]);
        }
        vector<mint> res(deg);
        res[0] = (*this)[0].inv();
        for (int i = 1; i < deg; i++) {
            mint r = 0;
            for (auto &&[k, val] : dat) {
                if (k > i) break;
                r -= val * res[i - k];
            }
            res[i] = r * res[0];
        }
        return res;
    }
    
    // log(f) = \int f'/f dx, f[0] must be 1
    constexpr FPS log(int deg = -1) const {
        assert(this->size() >= 1 && (*this)[0] == 1);
        if (count_terms() <= SPARSE_BOARDER) return log_sparse(deg);
        if (deg < 0) deg = (int)this->size();
        return ((diff() * inv(deg)).pre(deg - 1)).integral();
    }
    constexpr FPS log_sparse(int deg = -1) const {
        assert(this->size() >= 1 && (*this)[0] == 1);
        if (deg < 0) deg = (int)this->size();
        vector<pair<int, mint>> dat;
        for (int i = 1; i < (int)this->size(); i++) if ((*this)[i] != mint(0)) {
            dat.emplace_back(i, (*this)[i]);
        }
        BiCoef<mint> bc(deg);
        vector<mint> res(deg), tmp(deg);
        for (int i = 0; i < deg - 1; i++) {
            mint r = mint(i + 1) * (*this)[i + 1];
            for (auto &&[k, val] : dat) {
                if (k > i) break;
                r -= val * tmp[i - k];
            }
            tmp[i] = r;
            res[i + 1] = r * bc.inv(i + 1);
        }
        return res;
    }
    
    // exp(f), f[0] must be 0
    constexpr FPS exp(int deg = -1) const {
        if ((int)this->size() == 0) return {mint(1)};
        if (count_terms() <= SPARSE_BOARDER) return exp_sparse(deg);
        if constexpr (std::is_same_v<mint, Fp<998244353>>) return exp_ntt_friendly(deg);
        assert((*this)[0] == 0);
        if (deg < 0) deg = (int)this->size();
        FPS res(1, 1);
        for (int d = 1; d < deg; d <<= 1) {
            res = res * (pre(d << 1) - res.log(d << 1) + 1).pre(d << 1);
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS exp_ntt_friendly(int deg = -1) const {
        if ((int)this->size() == 0) return {mint(1)};
        assert((*this)[0] == 0);
        if (deg < 0) deg = (int)this->size();

        FPS fiv;
        fiv.reserve(deg + 1);
        fiv.emplace_back(mint(0));
        fiv.emplace_back(mint(1));

        auto inplace_integral = [&](FPS &F) -> void {
            const int n = (int)F.size();
            auto mod = mint::get_mod();
            while ((int)fiv.size() <= n) {
                int i = fiv.size();
                fiv.emplace_back((-fiv[mod % i]) * (mod / i));
            }
            F.insert(begin(F), mint(0));
            for (int i = 1; i <= n; i++) F[i] *= fiv[i];
        };

        auto inplace_diff = [](FPS &F) -> void {
            if (F.empty()) return;
            F.erase(begin(F));
            mint coef = 1;
            for (int i = 0; i < (int)F.size(); i++) {
                F[i] *= coef;
                coef++;
            }
        };

        FPS b{1, (1 < (int)this->size() ? (*this)[1] : 0)}, c{1}, z1, z2{1, 1};
        for (int m = 2; m < deg; m <<= 1) {
            auto y = b;
            y.resize(m * 2);
            ntt_trans(y);
            z1 = z2;
            FPS z(m);
            for (int i = 0; i < m; i++) z[i] = y[i] * z1[i];
            ntt_trans_inv(z);
            fill(begin(z), begin(z) + m / 2, mint(0));
            ntt_trans(z);
            for (int i = 0; i < m; i++) z[i] *= -z1[i];
            ntt_trans_inv(z);
            c.insert(end(c), begin(z) + m / 2, end(z));
            z2 = c;
            z2.resize(m * 2);
            ntt_trans(z2);
            FPS x(begin(*this), begin(*this) + min((int)this->size(), m));
            inplace_diff(x);
            x.emplace_back(mint(0));
            ntt_trans(x);
            for (int i = 0; i < m; i++) x[i] *= y[i];
            ntt_trans_inv(x);
            x -= b.diff();
            x.resize(m * 2);
            for (int i = 0; i < m - 1; i++) x[m + i] = x[i], x[i] = mint(0);
            ntt_trans(x);
            for (int i = 0; i < m * 2; i++) x[i] *= z2[i];
            ntt_trans_inv(x);
            x.pop_back();
            inplace_integral(x);
            for (int i = m; i < min((int)this->size(), m * 2); i++) x[i] += (*this)[i];
            fill(begin(x), begin(x) + m, mint(0));
            ntt_trans(x);
            for (int i = 0; i < m * 2; i++) x[i] *= y[i];
            ntt_trans_inv(x);
            b.insert(end(b), begin(x) + m, end(x));
        }
        return FPS(begin(b), begin(b) + deg);
    }
    constexpr FPS exp_sparse(int deg = -1) const {
        if ((int)this->size() == 0) return {mint(1)};
        assert((*this)[0] == 0);
        if (deg < 0) deg = (int)this->size();
        vector<pair<int, mint>> dat;
        for (int i = 1; i < (int)this->size(); i++) if ((*this)[i] != mint(0)) {
            dat.emplace_back(i - 1, (*this)[i] * i);
        }
        BiCoef<mint> bc(deg);
        vector<mint> res(deg);
        res[0] = 1;
        for (int i = 1; i < deg; i++) {
            mint r = 0;
            for (auto &&[k, val] : dat) {
                if (k > i - 1) break;
                r += val * res[i - k - 1];
            }
            res[i] = r * bc.inv(i);
        }
        return res;
    }
    
    // pow(f) = exp(e * log f)
    constexpr FPS pow(long long e, int deg = -1) const {
        if (count_terms() <= SPARSE_BOARDER) return pow_sparse(e, deg);
        assert(e >= 0);
        if (deg < 0) deg = (int)this->size();
        if (deg == 0) return FPS();
        if (e == 0) {
            FPS res(deg, 0);
            res[0] = 1;
            return res;
        }
        long long ord = 0;
        while (ord < (int)this->size() && (*this)[ord] == 0) ord++;
        if (ord == (int)this->size() || ord > (deg - 1) / e) return FPS(deg, 0);
        mint k = (*this)[ord];
        FPS res = ((((*this) >> ord) / k).log(deg) * e).exp(deg) * mint(k).pow(e) << (e * ord);
        res.resize(deg);
        return res;
    }
    constexpr FPS pow_sparse(long long e, int deg = -1) const {
        assert(e >= 0);
        if (deg < 0) deg = (int)this->size();
        if (deg == 0) return FPS();
        if (e == 0) {
            FPS res(deg, 0);
            res[0] = 1;
            return res;
        }
        long long ord = 0;
        while (ord < (int)this->size() && (*this)[ord] == 0) ord++;
        if (ord == (int)this->size() || ord > (deg - 1) / e) return FPS(deg, 0);
        if ((*this)[0] == 1) return pow_sparse_constant1(e, deg);
        auto f = (*this);
        rotate(f.begin(), f.begin() + ord, f.end());
        mint con = f[0], icon = f[0].inv();
        for (int i = 0; i < deg; i++) f[i] *= icon;
        auto res = f.pow_sparse_constant1(e, deg);
        int ord2 = e * ord;
        rotate(res.begin(), res.begin() + (deg - ord2), res.end());
        fill(res.begin(), res.begin() + ord2, mint(0));
        mint pw = con.pow(e);
        for (int i = ord2; i < deg; i++) res[i] *= pw;
        return res;
    }
    constexpr FPS pow_sparse_constant1(mint e, int deg = -1) const {
        assert((int)this->size() > 0 && (*this)[0] == 1);
        if (deg < 0) deg = (int)this->size();
        vector<pair<int, mint>> dat;
        for (int i = 1; i < (int)this->size(); i++) if ((*this)[i] != mint(0)) {
            dat.emplace_back(i, (*this)[i]);
        }
        BiCoef<mint> bc(deg);
        vector<mint> res(deg);
        res[0] = 1;
        for (int i = 0; i < deg - 1; i++) {
            mint &r = res[i + 1];
            for (auto &&[k, val] : dat) {
                if (k > i + 1) break;
                mint t = val * res[i - k + 1];
                r += t * (mint(k) * e - mint(i - k + 1));
            }
            r *= bc.inv(i + 1);
        }
        return res;
    }
    
    // sqrt(f)
    constexpr FPS sqrt(int deg = -1) const {
        if (count_terms() <= SPARSE_BOARDER) return sqrt_sparse(deg);
        if (deg < 0) deg = (int)this->size();
        if ((int)this->size() == 0) return FPS(deg, 0);
        if ((*this)[0] == mint(0)) {
            for (int i = 1; i < (int)this->size(); i++) {
                if ((*this)[i] != mint(0)) {
                    if (i & 1) return FPS();
                    if (deg - i / 2 <= 0) return FPS(deg, 0);
                    auto res = ((*this) >> i).sqrt(deg - i / 2);
                    if (res.empty()) return FPS();
                    res = res << (i / 2);
                    if ((int)res.size() < deg) res.resize(deg, mint(0));
                    return res;
                }
            }
            return FPS(deg, 0);
        }
        long long sqr = mod_sqrt<long long>((*this)[0].val, mint::get_mod());
        if (sqr == -1) return FPS();
        assert((*this)[0].val == sqr * sqr % mint::get_mod());
        FPS res = {mint(sqr)};
        mint iv2 = mint(2).inv();
        for (int d = 1; d < deg; d <<= 1) {
            res = (res + pre(d << 1) * res.inv(d << 1)).pre(d << 1) * iv2;
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS sqrt_sparse(int deg) const {
        if (deg < 0) deg = (int)this->size();
        if ((int)this->size() == 0) return FPS(deg, 0);
        if ((*this)[0] == mint(0)) {
            for (int i = 1; i < (int)this->size(); i++) {
                if ((*this)[i] != mint(0)) {
                    if (i & 1) return FPS();
                    if (deg - i / 2 <= 0) return FPS(deg, 0);
                    auto res = ((*this) >> i).sqrt_sparse(deg - i / 2);
                    if (res.empty()) return FPS();
                    res = res << (i / 2);
                    if ((int)res.size() < deg) res.resize(deg, mint(0));
                    return res;
                }
            }
            return FPS(deg, 0);
        }
        mint con = (*this)[0], icon = con.inv();
        long long sqr = mod_sqrt<long long>(con.val, mint::get_mod());
        if (sqr == -1) return FPS();
        assert(con.val == sqr * sqr % mint::get_mod());
        auto res = (*this) * icon;
        return res.sqrt_sparse_constant1(deg) * sqr;
    }
    constexpr FPS sqrt_sparse_constant1(int deg) const {
        return pow_sparse_constant1(mint(2).inv(), deg);
    }
    
    // friend operators
    friend constexpr FPS diff(const FPS &f) { return f.diff(); }
    friend constexpr FPS integral(const FPS &f) { return f.integral(); }
    friend constexpr FPS inv(const FPS &f, int deg = -1) { return f.inv(deg); }
    friend constexpr FPS log(const FPS &f, int deg = -1) { return f.log(deg); }
    friend constexpr FPS exp(const FPS &f, int deg = -1) { return f.exp(deg); }
    friend constexpr FPS pow(const FPS &f, long long e, int deg = -1) { return f.pow(e, deg); }
    friend constexpr FPS sqrt(const FPS &f, int deg = -1) { return f.sqrt(deg); }
};

// Bostan-Mori
// find [x^N] P(x)/Q(x), O(K log K log N)
// deg(Q(x)) = K, deg(P(x)) < K
template<typename mint> mint BostanMori
(const FPS<mint> &P, const FPS<mint> &Q, long long N) {
    assert(!P.empty() && !Q.empty());
    if (N == 0 || Q.size() == 1) return P[0] / Q[0];
    
    int qdeg = (int)Q.size();
    FPS<mint> P2{P}, minusQ{Q};
    P2.resize(qdeg - 1);
    for (int i = 1; i < (int)Q.size(); i += 2) minusQ[i] = -minusQ[i];
    P2 *= minusQ;
    FPS<mint> Q2 = Q * minusQ;
    FPS<mint> S(qdeg - 1), T(qdeg);
    for (int i = 0; i < (int)S.size(); ++i) {
        S[i] = (N % 2 == 0 ? P2[i * 2] : P2[i * 2 + 1]);
    }
    for (int i = 0; i < (int)T.size(); ++i) {
        T[i] = Q2[i * 2];
    }
    return BostanMori(S, T, N >> 1);
}

// find generating function F(x) (= P(x) / Q(x)) of linearly D-recurrent sequence, O(D log D)
// x[0] = A[0], x[1] = A[1], ..., x[D-1] = A[D-1]
// x[i] = C[0]x[i-1] + C[1]x[i-2] + ... + C[D-1]x[i-D]
template<typename mint> pair<FPS<mint>, FPS<mint>> find_generating_function
(const vector<mint> &A, const vector<mint> &C) {
    assert(A.size() == C.size());
    int D = (int)C.size();
    FPS<mint> Q(D+1);
    Q[0] = 1;
    for (int i = 1; i <= D; i++) Q[i] = -C[i-1];
    FPS<mint> P = (Q * FPS<mint>(A)).pre(D);
    return make_pair(P, Q);  
}

// find x[K] of linearly D-recurrent sequence, O(D log D log K)
// x[0] = A[0], x[1] = A[1], ..., x[D-1] = A[D-1]
// x[i] = C[0]x[i-1] + C[1]x[i-2] + ... + C[D-1]x[i-D]
template<typename mint> mint kth_term
(const vector<mint> &A, const vector<mint> &C, long long K) {
    assert(A.size() == C.size());
    auto [P, Q] = find_generating_function(A, C);
    return BostanMori(P, Q, K);  
}

// Berlekamp-Massey, find linear recurrence, O(D^2)
// given: A[0], ..., A[N-1]
// find: C[0], ..., C[D-1] s.t. A[i] = C[0]A[i-1] + C[1]x[i-2] + ... + C[D-1]x[i-D]
template<typename mint> vector<mint> BerlekampMassey(const vector<mint> &A) {
    const int N = (int)A.size();
    vector<mint> b({mint(-1)}), c({mint(-1)});
    mint x = 0, y = 1;
    for (int iter = 1; iter <= N; iter++) {
        int cl = c.size(), bl = b.size();
        x = 0;
        for (int i = 0; i < cl; i++) x += c[i] * A[iter - cl + i];
        b.emplace_back(0);
        bl++;
        if (x == mint(0)) continue;
        mint f = x / y;
        if (cl < bl) {
            auto tmp = c;
            c.insert(c.begin(), bl - cl, mint(0));
            for (int i = 0; i < bl; i++) c[bl - i - 1] -= f * b[bl - i - 1];
            b = tmp;
            y = x;
        } else {
            for (int i = 0; i < bl; i++) c[cl - i - 1] -= f * b[bl - i - 1];
        }
    }
    c.pop_back();
    reverse(c.begin(), c.end());
    return c;
}
template<typename mint> pair<FPS<mint>, FPS<mint>> find_generating_function
(const vector<mint> &A) {
    auto C = BerlekampMassey(A);
    vector<mint> A2;
    for (int i = 0; i < (int)C.size(); i++) A2.emplace_back(A[i]);
    return find_generating_function(A2, C);
}
template<typename mint> mint BMBM(const vector<mint> &A, long long K) {
    auto [P, Q] = find_generating_function(A);
    return BostanMori(P, Q, K);   
}

// composition of FPS, calc g(f(x)), O(N (log N)^2)
template<class mint>
FPS<mint> composition(FPS<mint> g, FPS<mint> f, int deg = -1) {
    auto rec = [&](auto &&rec, FPS<mint> Q, int n, int h, int k) -> FPS<mint> {
        if (n == 0) {
            FPS<mint> T{begin(Q), begin(Q) + k};
            T.emplace_back(mint(1));
            FPS<mint> u = g * T.rev().inv().rev();
            FPS<mint> P(h * k);
            for (int i = 0; i < (int)g.size(); i++) P[k - i - 1] = u[i + k];
            return P;
        }
        FPS<mint> nQ(h * k * 4), nR(h * k * 2);
        for (int i = 0; i < k; i++) {
            copy(begin(Q) + i * h, begin(Q) + i * h + n + 1, begin(nQ) + i * h * 2);
        }
        nQ[h * k * 2] += 1;
        ntt_trans(nQ);
        for (int i = 0; i < h * k * 4; i += 2) swap(nQ[i], nQ[i + 1]);
        for (int i = 0; i < h * k * 2; i++) nR[i] = nQ[i * 2] * nQ[i * 2 + 1];
        ntt_trans_inv(nR);
        nR[0] -= 1;
        Q.assign(h * k, 0);
        for (int i = 0; i < k * 2; i++) for (int j = 0; j <= n / 2; j++) {
            Q[i * h / 2 + j] = nR[i * h + j];
        }
        auto P = rec(rec, Q, n / 2, h / 2, k * 2);
        FPS<mint> nP(h * k * 4);
        for (int i = 0; i < k * 2; i++) for (int j = 0; j <= n / 2; j++) {
            nP[i * h * 2 + j * 2 + n % 2] = P[i * h / 2 + j];
        }
        ntt_trans(nP);
        for (int i = 1; i < h * k * 4; i <<= 1) reverse(begin(nQ) + i, begin(nQ) + i * 2);
        for (int i = 0; i < h * k * 4; i++) nP[i] *= nQ[i];
        ntt_trans_inv(nP);
        P.assign(h * k, 0);
        for (int i = 0; i < k; i++) {
            copy(begin(nP) + i * h * 2, begin(nP) + i * h * 2 + n + 1, begin(P) + i * h);
        }
        return P;
    };
    if (deg == -1) deg = max((int)f.size(), (int)g.size());
    f.resize(deg), g.resize(deg);
    int n = (int)f.size() - 1, h = 1, k = 1;
    while (h < n + 1) h *= 2;
    FPS<mint> Q(h * k);
    for (int i = 0; i <= n; i++) Q[i] = -f[i];
    FPS<mint> P = rec(rec, Q, n, h, k);
    return P.pre(n + 1).rev();
}

// Power Projection, O(N (log N)^2)
// for i = 0, 1, ..., m, calc [x^(f の最高次数)] f(x)^i g(x) 
template<class mint, int MOD = mint::get_mod(), int pr = calc_primitive_root(MOD)>
FPS<mint> power_projection(FPS<mint> f, FPS<mint> g = {1}, int m = -1) {
    int n = (int)f.size() - 1, k = 1, h = 1;
    g.resize(n + 1);
    if (m < 0) m = n;
    while (h < n + 1) h <<= 1;
    FPS<mint> P((n + 1) * k), Q((n + 1) * k), nP, nQ, buf, buf2;
    for (int i = 0; i <= n; i++) P[i * k] = g[i];
    for (int i = 0; i <= n; i++) Q[i * k] = -f[i];
    Q[0]++;
    mint iv2 = mint(2).inv();
    while (n) {
        mint w = mint(pr).pow((MOD - 1) / (2 * k)), iw = w.inv();
        buf2.resize(k);
        auto ntt_doubling = [&]() {
            copy(begin(buf), end(buf), begin(buf2));
            ntt_trans_inv(buf2);
            mint c = 1;
            for (int i = 0; i < k; i++) buf2[i] *= c, c *= w;
            ntt_trans(buf2);
            copy(begin(buf2), end(buf2), back_inserter(buf));
        };
        nP.clear(), nQ.clear();
        for (int i = 0; i <= n; i++) {
            buf.resize(k);
            copy(begin(P) + i * k, begin(P) + (i + 1) * k, begin(buf));
            ntt_doubling();
            copy(begin(buf), end(buf), back_inserter(nP));
            buf.resize(k);
            copy(begin(Q) + i * k, begin(Q) + (i + 1) * k, begin(buf));
            if (i == 0) {
                for (int j = 0; j < k; j++) buf[j]--;
                ntt_doubling();
                for (int j = 0; j < k; j++) buf[j]++;
                for (int j = 0; j < k; j++) buf[k + j]--;
            } else {
                ntt_doubling();
            }
            copy(begin(buf), end(buf), back_inserter(nQ));
        }
        nP.resize(h * 2 * k * 2), nQ.resize(h * 2 * k * 2);
        FPS<mint> p(h * 2), q(h * 2);
        w = mint(pr).pow((MOD - 1) / (h * 2)), iw = w.inv();
        vector<int> btr;
        if (n % 2) {
            btr.resize(h);
            for (int i = 0, lg = bsf(h); i < h; i++) {
                btr[i] = (btr[i >> 1] >> 1) + ((i & 1) << (lg - 1));
            }
        }
        for (int j = 0; j < k * 2; j++) {
            p.assign(h * 2, 0), q.assign(h * 2, 0);
            for (int i = 0; i < h; i++) p[i] = nP[i * k * 2 + j], q[i] = nQ[i * k * 2 + j];
            ntt_trans(p), ntt_trans(q);
            for (int i = 0; i < h * 2; i += 2) swap(q[i], q[i + 1]);
            for (int i = 0; i < h * 2; i++) p[i] *= q[i];
            for (int i = 0; i < h; i++) q[i] = q[i * 2] * q[i * 2 + 1];
            if (n & 1) {
                mint c = iv2;
                buf.resize(h);
                for (int i : btr) buf[i] = (p[i * 2] - p[i * 2 + 1]) * c, c *= iw;
                swap(p, buf);
            } else {
                for (int i = 0; i < h; i++) p[i] = (p[i * 2] + p[i * 2 + 1]) * iv2;
            }
            p.resize(h), q.resize(h);
            ntt_trans_inv(p), ntt_trans_inv(q);
            for (int i = 0; i < h; i++) nP[i * k * 2 + j] = p[i];
            for (int i = 0; i < h; i++) nQ[i * k * 2 + j] = q[i];
        }
        nP.resize((n / 2 + 1) * k * 2), nQ.resize((n / 2 + 1) * k * 2);
        swap(P, nP), swap(Q, nQ);
        n /= 2, h /= 2, k *= 2;
    }
    FPS<mint> S{begin(P), begin(P) + k}, T{begin(Q), begin(Q) + k};
    ntt_trans_inv(S), ntt_trans_inv(T);
    T[0]--;
    if (T[0] == 0) return S.rev().pre(m + 1);
    else return (S.rev() * (T + (FPS<mint>{1} << k)).rev().inv(m + 1)).pre(m + 1);
}

// find g s.t. f(g(x)) ≡ x (mod x^{deg}), O(N (log N)^2)
template<class mint>
FPS<mint> compositional_inverse(FPS<mint> f, int deg = -1) {
    assert((int)f.size() >= 2 && f[1] != 0);
    if (deg == -1) deg = (int)f.size();
    if (deg < 2) return FPS<mint>{0, f[1].inv()}.pre(deg);
    int n = deg - 1;
    FPS<mint> h = power_projection(f) * n;
    for (int k = 1; k <= n; k++) h[k] /= k;
    h = h.rev(), h *= h[0].inv();
    FPS<mint> g = (h.log() * mint(-n).inv()).exp();
    g *= f[1].inv();
    return (g << 1).pre(deg);
}


//------------------------------//
// Polynomial Algorithms
//------------------------------//

// all inverse
template<class mint> vector<mint> all_inverse(const vector<mint> &v) {
    for (auto &&vi : v) assert(vi != mint(0));
    int N = (int)v.size();
    vector<mint> res(N + 1, mint(1));
    for (int i = 0; i < N; i++) res[i + 1] = res[i] * v[i];
    mint t = res.back().inv();
    res.pop_back();
    for (int i = N - 1; i >= 0; i--) res[i] *= t, t *= v[i];
    return res;
}

// polynomial merge technique
template<class mint> FPS<mint> all_product(const vector<FPS<mint>> &fs) {
    using Node = pair<int, FPS<mint>>;
    priority_queue<Node, vector<Node>, greater<Node>> que;
    que.push({0, FPS<mint>{1}});
    for (const auto &f : fs) que.push({(int)f.size(), f});
    while (que.size() >= 2) {
        auto [df, f] = que.top();
        que.pop();
        auto [dg, g] = que.top();
        que.pop();
        auto h = f * g;
        que.push({(int)h.size(), h});
    }
    auto [deg, res] = que.top();
    return res;
}

// find f(x)^n mod g(x)
template<class mint, class T_VAL = long long> 
FPS<mint> mod_pow(const FPS<mint> &f, T_VAL e, const FPS<mint> &mod) {
    assert(!mod.empty());
    auto iv = mod.rev().inv();
    auto calc_quo = [&](const FPS<mint> &pol) -> FPS<mint> {
        if (pol.size() < mod.size()) return FPS<mint>();
        int deg = (int)pol.size() - (int)mod.size() + 1;
        return (pol.rev().pre(deg) * iv.pre(deg)).pre(deg).rev();
    };
    FPS<mint> res{1}, b(f);
    while (e) {
        if (e & 1) res *= b, res -= calc_quo(res) * mod;
        b *= b;
        b -= calc_quo(b) * mod;
        e >>= 1;
        assert(b.size() + 1 <= mod.size());
        assert(res.size() + 1 <= mod.size());
    }
    return res;
}

// middle product 
// c[i] = sum_j a[i+j]b[j] (deg(a) = n, deg(b) = m -> deg(c) = n-m)
template<class mint>
FPS<mint> middle_product(const FPS<mint> &a, const FPS<mint> &b) {
    assert(a.size() >= b.size());
    if (b.empty()) return FPS<mint>((int)a.size() - (int)b.size() + 1);
    int N = 1;
    while (N < (int)a.size()) N <<= 1;
    FPS<mint> fa(N), fb(N);
    copy(a.begin(), a.end(), fa.begin());
    copy(b.rbegin(), b.rend(), fb.begin());
    fa *= fb;
    fa.resize(a.size());
    fa.erase(fa.begin(), fa.begin() + (int)b.size() - 1);
    return fa;
}

// Polynomial Taylor Shift
template<class mint> FPS<mint> taylor_shift(const FPS<mint> &f, long long c) {
    int N = (int)f.size() - 1;
    BiCoef<mint> bc(N + 1);
    FPS<mint> p(N + 1), q(N + 1);
    for (int i = 0; i <= N; i++) {
        p[i] = f[i] * bc.fact(i);
        q[N - i] = mint(c).pow(i) * bc.finv(i);
    }
    FPS<mint> pq = p * q;
    FPS<mint> res(N + 1);
    for (int i = 0; i <= N; i++) res[i] = pq[i + N] * bc.finv(i);
    return res;
}

// multipoint evaluation, polynomial interpolation
template<class mint> struct SubproductTree {
    // inner data
    int num_points, siz;
    vector<FPS<mint>> tree;

    // constructor
    SubproductTree() {}
    SubproductTree(const vector<mint> &x) {
        num_points = (int)x.size();
        siz = 1;
        while (siz < num_points) siz <<= 1;
        tree.resize(siz * 2);
        for (int i = 0; i < siz; i++) tree[siz + i] = {1, (i < num_points ? -x[i] : 0)};
        for (int i = siz - 1; i >= 1; i--) tree[i] = tree[i * 2] * tree[i * 2 + 1];
    }

    // multipoint evaluation
    vector<mint> eval(FPS<mint> f) {
        int N = (int)f.size();
        if (N == 0) return vector<mint>(num_points, mint(0));
        f.resize(N * 2 - 1);
        vector<FPS<mint>> g(siz * 2);
        g[1] = tree[1];
        g[1].resize(N);
        g[1] = inv(g[1]);
        g[1] = middle_product(f, g[1]);
        g[1].resize(siz);
        for (int i = 1; i < siz; i++) {
            g[i * 2] = middle_product(g[i], tree[i * 2 + 1]);
            g[i * 2 + 1] = middle_product(g[i], tree[i * 2]);
        }
        vector<mint> res(num_points);
        for (int i = 0; i < num_points; i++) res[i] = g[siz + i][0];
        return res;
    }

    // polynomial interpolation
    FPS<mint> interpolate(const vector<mint> &y) {
        assert((int)y.size() == num_points);
        vector<mint> p(num_points);
        for (int i = 0; i < num_points; i++) p[i] = tree[1][num_points - i - 1] * (i + 1);
        p = eval(p);
        vector<FPS<mint>> t(siz * 2);
        for (int i = 0; i < siz; i++) t[siz + i] = {(i < num_points ? y[i] / p[i] : 0)};
        for (int i = siz - 1; i >= 1; i--) {
            t[i] = t[i * 2] * tree[i * 2 + 1];
            auto rt = t[i * 2 + 1] * tree[i * 2];
            for (int k = 0; k < (int)t[i].size(); k++) t[i][k] += rt[k];
        }
        t[1].resize(num_points);
        reverse(t[1].begin(), t[1].end());
        return t[1];
    }
};

// multipoint evaluation, polynomial interpolation
template<class mint>
vector<mint> multipoint_eval(const FPS<mint> &f, const vector<mint> &x) {
    if (x.empty()) return {};
    SubproductTree<mint> st(x);
    return st.eval(f);
}
template<class mint>
FPS<mint> interpolate(const vector<mint> &x, const vector<mint> &y) {
    assert(x.size() == y.size());
    if (x.empty()) return {};
    SubproductTree<mint> st(x);
    return st.interpolate(y);
}

// multipoint evaluation (case: geometric sequence)
// for k = 0, 1, ..., M-1, calc f(ar^k)
template<class mint>
vector<mint> multipoint_eval(const FPS<mint> &f, const mint &a, const mint &r, int M) {
    // calc 1, 1, r, r^3, r^6, r^10, ...
    auto calc = [&](const mint &r, int m) -> FPS<mint> {
        FPS<mint> res(m, mint(1));
        mint po = 1;
        for (int i = 0; i < m - 1; i++) res[i + 1] = res[i] * po, po *= r;
        return res;
    };
    int N = (int)f.size();
    if (M == 0) return vector<mint>();
    if (r == mint(0)) {
        vector<mint> res(M);
        for (int i = 1; i < M; i++) res[i] = f[0];
        res[0] = f.eval(a);
        return res;
    }
    if (min(N, M) < 60) {
        vector<mint> res(M);
        mint b = a;
        for (int i = 0; i < M; i++) res[i] = f.eval(b), b *= r;
        return res;
    }
    FPS<mint> res = f;
    mint po = 1;
    for (int i = 0; i < N; i++) res[i] *= po, po *= a;
    FPS<mint> A = calc(r, N + M - 1), B = calc(r.inv(), max(N, M));
    for (int i = 0; i < N; i++) res[i] *= B[i];
    res = middle_product(A, res);
    for (int i = 0; i < M; i++) res[i] *= B[i];
    return res;
}

// polynomial interpolation (case: geometric sequence)
// y[i] = f(ar^i) -> find f
template<class mint>
vector<mint> interpolate(const mint &a, const mint &r, const FPS<mint> &y) {
    int N = (int)y.size();
    if (N == 0) return FPS<mint>();
    if (N == 1) return {y[0]};
    auto Y = y;
    mint ir = r.inv(), ia = a.inv();
    FPS<mint> po(N + N - 1, 1), po2(N + N - 1, 1), ipo(N + N - 1, 1), ipo2(N + N - 1, 1);
    for (int i = 0; i < N + N - 2; i++) po[i + 1] = po[i] * r, po2[i + 1] = po2[i] * po[i];
    for (int i = 0; i < N; i++) ipo[i + 1] = ipo[i] * ir, ipo2[i + 1] = ipo2[i] * ipo[i];
    vector<mint> S(N, mint(1));
    for (int i = 1; i < N; i++) S[i] = S[i - 1] * (mint(1) - po[i]);
    vector<mint> iS = all_inverse(S);
    mint sn = S[N - 1] * (mint(1) - po[N]);
    for (int i = 0; i < N; i++) {
        Y[i] = Y[i] * po2[N - i - 1] * ipo2[N - 1] * iS[i] * iS[N - i - 1];
        if (i & 1) Y[i] = -Y[i];
    }
    for (int i = 0; i < N; i++) Y[i] *= ipo2[i];
    FPS<mint> f = middle_product(po2, Y);
    for (int i = 0; i < N; i++) f[i] *= ipo2[i];
    FPS<mint> g(N, mint(1));
    for (int i = 1; i < N; i++) {
        g[i] = po2[i] * sn * iS[i] * iS[N - i];
        if (i & 1) g[i] = -g[i];
    }
    f = f * g;
    f.resize(N);
    reverse(f.begin(), f.end());
    mint p = 1;
    for (int i = 0; i < N; i++) f[i] *= p, p *= ia;
    return f;
}

// shift of sampling points
// input y[i] = f(i) (i = 0, 1, ..., N-1)
// output f[c], f[c+1], ..., f[c+M-1]
template<class mint> vector<mint> shift_of_sampling_points(const vector<mint> &y, const mint &c, int M = -1) {
    int N = (int)y.size();
    if (M == -1) M = N;
    int T = c.val;
    if (T < N) {
        FPS<mint> res(M);
        int ptr = 0;
        for (int i = T; i < N && ptr < M; i++) res[ptr++] = y[i];
        if (N < T + M) {
            auto suffix = shift_of_sampling_points(y, mint(N), M - ptr);
            for (int i = N; i < T + M; i++) res[ptr++] = suffix[i - N];
        }
        return res;
    }
    if (T + M > mint::get_mod()) {
        auto prefix = shift_of_sampling_points(y, mint(T), mint::get_mod() - T);
        auto suffix = shift_of_sampling_points(y, mint(0), M - (int)prefix.size());
        copy(begin(suffix), end(suffix), back_inserter(prefix));
        return prefix;
    }
    BiCoef<mint> bc(N);
    FPS<mint> res(M), d(N), h(M + N - 1);
    for (int i = 0; i < N; i++) {
        d[i] = bc.finv(i) * bc.finv(N - i - 1) * y[i];
        if ((N - i - 1) & 1) d[i] = -d[i];
    }
    for (int i = 0; i < M + N - 1; i++) h[i] = mint(T - N + i + 1);
    h = all_inverse(h);
    auto dh = d * h;
    mint cur = T;
    for (int i = 1; i < N; i++) cur *= T - i;
    for (int i = 0; i < M; i++) {
        res[i] = cur * dh[N + i - 1];
        cur *= T + i + 1;
        cur *= h[i];
    }
    return res;
}


//------------------------------//
// Matrix
//------------------------------//

// modint matrix
template<class mint> struct MintMatrix {
    // inner value
    int H, W;
    vector<vector<mint>> val;
    
    // constructors
    MintMatrix() {}
    MintMatrix(const MintMatrix&) = default;
    MintMatrix& operator = (const MintMatrix&) = default;
    MintMatrix(int h, int w) : H(h), W(w), val(h, vector<mint>(w)) {}
    MintMatrix(int h, int w, mint x) : H(h), W(w), val(h, vector<mint>(w, x)) {}
    void init(int h, int w, mint x) {
        H = h, W = w;
        val.assign(h, vector<mint>(w, x));
    }
    void resize(int h, int w) {
        H = h, W = w;
        val.resize(h);
        for (int i = 0; i < h; ++i) val[i].resize(w);
    }
    
    // getter and debugger
    constexpr int height() const { return H; }
    constexpr int width() const { return W; }
    constexpr bool empty() const { return height() == 0; }
    vector<mint>& operator [] (int i) { return val[i]; }
    const vector<mint>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const MintMatrix &mat) {
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) {
                if (j) os << ' ';
                os << mat.val[i][j];
            }
            os << '\n';
        }
        return os;
    }
    
    // comparison operators
    constexpr bool operator == (const MintMatrix &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const MintMatrix &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr MintMatrix& operator += (const MintMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] + r.val[i][j];
        return *this;
    }
    constexpr MintMatrix& operator -= (const MintMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] - r.val[i][j];
        return *this;
    }
    constexpr MintMatrix& operator *= (const mint &v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] * v;
        return *this;
    }
    constexpr MintMatrix& operator *= (const MintMatrix &r) {
        assert(width() == r.height());
        MintMatrix<mint> res(height(), r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] = res[i][j] + val[i][k] * r.val[k][j];
        return (*this) = res;
    }
    constexpr MintMatrix operator + () const { 
        return MintMatrix(*this);
    }
    constexpr MintMatrix operator - () const {
        MintMatrix res(*this);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res.val[i][j] = -res.val[i][j];
        return res;
    }
    constexpr MintMatrix operator + (const MintMatrix &r) const { 
        return MintMatrix(*this) += r;
    }
    constexpr MintMatrix operator - (const MintMatrix &r) const {
        return MintMatrix(*this) -= r;
    }
    constexpr MintMatrix operator * (const mint &v) const {
        return MintMatrix(*this) *= v;
    }
    constexpr MintMatrix operator * (const MintMatrix &r) const {
        return MintMatrix(*this) *= r;
    }
    constexpr vector<mint> operator * (const vector<mint> &v) const {
        assert(width() == v.size());
        vector<mint> res(height(), mint(0));
        for (int i = 0; i < height(); i++)
            for (int j = 0; j < width(); j++)
                res[i] += val[i][j] * v[j];
        return res;
    }

    // transpose
    constexpr MintMatrix trans() const {
        MintMatrix<mint> res(width(), height());
        for (int row = 0; row < width(); row++)
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend constexpr MintMatrix trans(const MintMatrix &mat) {
        return mat.trans();
    }
    
    // pow
    constexpr MintMatrix pow(long long n) const {
        assert(height() == width());
        MintMatrix<mint> res(height(), width());
        MintMatrix<mint> mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = mint(1);
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr MintMatrix pow(const MintMatrix &mat, long long n) {
        return mat.pow(n);
    }
    
    // gauss-jordan
    constexpr int find_pivot(int cur_rank, int col) const {
        int pivot = -1;
        for (int row = cur_rank; row < height(); ++row) {
            if (val[row][col] != mint(0)) {
                pivot = row;
                break;
            }
        }
        return pivot;
    }
    constexpr void sweep(int cur_rank, int col, int pivot, bool sweep_upper = true) {
        swap(val[pivot], val[cur_rank]);
        auto ifac = val[cur_rank][col].inv();
        for (int col2 = cur_rank; col2 < width(); ++col2) {
            val[cur_rank][col2] *= ifac;
        }
        int row_start = (sweep_upper ? 0 : cur_rank + 1);
        for (int row = row_start; row < height(); ++row) {
            if (row != cur_rank && val[row][col] != mint(0)) {
                auto fac = val[row][col];
                for (int col2 = cur_rank; col2 < width(); ++col2) {
                    val[row][col2] -= val[cur_rank][col2] * fac;
                }
            }
        }
    }
    constexpr int gauss_jordan(int not_sweep_width = 0, bool sweep_upper = true) {
        int rank = 0;
        for (int col = 0; col < width(); ++col) {
            if (col == width() - not_sweep_width) break;
            int pivot = find_pivot(rank, col);
            if (pivot == -1) continue;
            sweep(rank++, col, pivot, sweep_upper);
        }
        return rank;
    }
    constexpr int gauss_jordan(vector<int> &core, int not_sweep_width, bool sweep_upper = true) {
        core.clear();
        int rank = 0;
        for (int col = 0; col < width(); ++col) {
            if (col == width() - not_sweep_width) break;
            int pivot = find_pivot(rank, col);
            if (pivot == -1) continue;
            core.push_back(col);
            sweep(rank++, col, pivot, sweep_upper);
        }
        return rank;
    }
    friend constexpr int gauss_jordan(MintMatrix &mat, int not_sweep_width = 0, bool sweep_upper = true) {
        return mat.gauss_jordan(not_sweep_width, sweep_upper);
    }

    // rank
    constexpr int get_rank() const {
        if (height() == 0 || width() == 0) return 0;
        MintMatrix A(*this);
        if (height() < width()) A = A.trans();
        return A.gauss_jordan(0, false);
    }
    friend constexpr int get_rank(const MintMatrix &mat) {
        return mat.get_rank();
    }

    // find one solution
    friend constexpr int linear_equation
    (const MintMatrix &mat, const vector<mint> &b, vector<mint> &res) {
        // extend
        MintMatrix<mint> A(mat.height(), mat.width() + 1);
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) A[i][j] = mat.val[i][j];
            A[i].back() = b[i];
        }
        int rank = A.gauss_jordan(1);
        
        // check if it has no solution
        for (int row = rank; row < mat.height(); ++row) if (A[row].back() != 0) return -1;

        // answer
        res.assign(mat.width(), 0);
        for (int i = 0; i < rank; ++i) res[i] = A[i].back();
        return rank;
    }
    friend constexpr int linear_equation(const MintMatrix &mat, const vector<mint> &b) {
        vector<mint> res;
        return linear_equation(mat, b, res);
    }

    // find all solutions
    friend int linear_equation
    (const MintMatrix &mat, const vector<mint> &b, vector<mint> &res, vector<vector<mint>> &zeros) {
        // extend
        MintMatrix<mint> A(mat.height(), mat.width() + 1);
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) A[i][j] = mat.val[i][j];
            A[i].back() = b[i];
        }
        vector<int> core;
        int rank = A.gauss_jordan(core, 1);
        
        // check if it has no solution
        for (int row = rank; row < mat.height(); ++row) {
            if (A[row].back() != 0) return -1;
        }

        // construct the core solution
        res.assign(mat.width(), mint(0));
        for (int i = 0; i < (int)core.size(); i++) res[core[i]] = A[i].back();
    
        // construct the all solutions
        zeros.clear();
        vector<bool> use(mat.width(), 0);
        for (auto c : core) use[c] = true;
        for (int j = 0; j < mat.width(); j++) {
            if (use[j]) continue;
            vector<mint> zero(mat.width(), mint(0));
            zero[j] = mint(1);
            for (int i = 0; i < (int)core.size(); i++) zero[core[i]] = -A[i][j];
            zeros.push_back(zero);
        }
        return rank;
    }
    
    // determinant
    constexpr mint det() const {
        assert(height() == width());
        if (height() == 0) return mint(1);
        MintMatrix<mint> A(*this);
        int rank = 0;
        mint res = mint(1);
        for (int col = 0; col < width(); ++col) {
            int pivot = A.find_pivot(rank, col);
            if (pivot == -1) return mint(0);
            if (pivot != rank) res = -res;
            res *= A[pivot][rank];
            A.sweep(rank++, col, pivot, false);
        }
        return res;
    }
    friend constexpr mint det(const MintMatrix &mat) {
        return mat.det();
    }
    constexpr mint det_nonprime_mod() const {
        assert(height() == width());
        if (height() == 0) return mint(1);
        MintMatrix<mint> A(*this);
        int rank = 0;
        mint res = mint(1);
        for (int col = 0; col < width(); ++col) {
            int pivot = A.find_pivot(rank, col);
            if (pivot == -1) return mint(0);
            if (pivot != rank) swap(A[pivot], A[rank]), res = -res;
            for (int row = rank + 1; row < height(); ++row) {
                while (A[row][col] != 0) {
                    swap(A[rank], A[row]), res = -res;
                    long long quo = A[row][col].get() / A[rank][col].get();
                    for (int col2 = rank; col2 < width(); ++col2) {
                        A[row][col2] -= A[rank][col2] * quo;
                    }
                }
            }
            rank++;
        }
        for (int col = 0; col < height(); ++col) res *= A[col][col];
        return res;
    }
    friend constexpr mint det_nonprime_mod(const MintMatrix &mat) {
        return mat.det_nonprime_mod();
    }

    // inv
    constexpr MintMatrix inv() const {
        assert(height() == width());

        // extend
        MintMatrix<mint> A(height(), width() + height());
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) A[i][j] = val[i][j];
            A[i][i+width()] = mint(1);
        }
        vector<int> core;
        int rank = A.gauss_jordan(height(), true);

        // gauss jordan
        if (rank < height()) return MintMatrix();
        MintMatrix<mint> res(height(), width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res[i][j] = A[i][j+width()];
        return res;
    }
    friend constexpr MintMatrix inv(const MintMatrix &mat) {
        return mat.inv();
    }
};


//------------------------------//
// Flow
//------------------------------//

// edge class (for max-flow)
template<class FLOW> struct FlowEdge {
    // core members
    int rev, from, to;
    FLOW cap, icap, flow;
    
    // constructor
    FlowEdge() {}
    FlowEdge(int r, int f, int t, FLOW c) : rev(r), from(f), to(t), cap(c), icap(c), flow(0) {}
    void reset() { cap = icap, flow = 0; }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowEdge& E) {
        return s << E.from << "->" << E.to << '(' << E.flow << '/' << E.icap << ')';
    }
};

// graph class (for max-flow)
template<class FLOW> struct FlowGraph {
    // core members
    vector<vector<FlowEdge<FLOW>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    
    // constructor
    FlowGraph(int n = 0) : list(n) { }
    void init(int n = 0) {
        list.clear(), list.resize(n);
        pos.clear();
    }
    void clear() {
        list.clear(), pos.clear();
    }
    
    // getter
    vector<FlowEdge<FLOW>> &operator [] (int i) {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    const vector<FlowEdge<FLOW>> &operator [] (int i) const {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowEdge<FLOW> &get_rev_edge(const FlowEdge<FLOW> &e) {
        return list[e.to][e.rev];
    }
    FlowEdge<FLOW> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowEdge<FLOW> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowEdge<FLOW>> get_edges() const {
        vector<FlowEdge<FLOW>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowEdge<FLOW> &e : list[i]) e.reset();
        }
    }
    void change_edge(FlowEdge<FLOW> &e, FLOW new_cap, FLOW new_flow) {
        assert(0 <= new_flow && new_flow <= new_cap);
        FlowEdge<FLOW> &re = get_rev_edge(e);
        e.cap = new_cap - new_flow, e.icap = new_cap, e.flow = new_flow;
        re.cap = new_flow;
    }
    
    // add_edge
    void add_edge(int from, int to, FLOW cap, FLOW rcap = 0) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowEdge<FLOW>(to_id, from, to, cap));
        list[to].push_back(FlowEdge<FLOW>(from_id, to, from, rcap));
    }

    // augment
    FLOW augment(int s, int t, FLOW up_flow = numeric_limits<FLOW>::max()) {
        vector<bool> seen(size(), false);
        auto dfs = [&](auto &&dfs, int v, FLOW up_flow) -> FLOW {
            if (v == t) return up_flow;
            seen[v] = true;
            for (int i = 0; i < (int)list[v].size(); i++) {
                FlowEdge<FLOW> &e = list[v][i], &re = get_rev_edge(e);
                if (seen[e.to] || e.cap <= 0) continue;
                FLOW flow = dfs(dfs, e.to, min(up_flow, e.cap));
                if (flow > 0) {
                    e.cap -= flow, e.flow += flow;
                    re.cap += flow, re.flow -= flow;
                    return flow;
                }
            }  
            return FLOW(0); 
        };
        return dfs(dfs, s, up_flow);
    };

    // find reachable nodes from node s (1: s-domain, -1: t-domain, 0: no reach)
    vector<int> find_cut(int s, int t) {
        vector<int> res(size(), 0);
        auto dfs_s = [&](auto &&dfs_s, int v) -> void {
            res[v] = 1;
            for (const auto &e : list[v]) {
                if (res[e.to] || e.cap <= 0) continue;
                dfs_s(dfs_s, e.to);
            }
        };
        auto dfs_t = [&](auto &&dfs_t, int v) -> void {
            res[v] = -1;
            for (const auto &e : list[v]) {
                auto re = get_rev_edge(e);
                if (res[e.to] || re.cap <= 0) continue;
                dfs_t(dfs_t, e.to);
            }
        };
        dfs_s(dfs_s, s), dfs_t(dfs_t, t);
        return res;
    }

    // debug
    friend ostream& operator << (ostream& s, const FlowGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// Dinic
template<class FLOW> FLOW Dinic(FlowGraph<FLOW> &G, int s, int t, FLOW limit_flow) {
    assert(0 <= s && s < G.size() && 0 <= t && t < G.size() && s != t);
    FLOW current_flow = 0;
    vector<int> level((int)G.size(), -1), iter((int)G.size(), 0);
    
    // Dinic BFS
    auto bfs = [&]() -> void {
        level.assign((int)G.size(), -1);
        level[s] = 0;
        queue<int> que;
        que.push(s);
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            for (const FlowEdge<FLOW> &e : G[v]) {
                if (level[e.to] < 0 && e.cap > 0) {
                    level[e.to] = level[v] + 1;
                    if (e.to == t) return;
                    que.push(e.to);
                }
            }
        }
    };
    
    // Dinic DFS
    auto dfs = [&](auto self, int v, FLOW up_flow) {
        if (v == t) return up_flow;
        FLOW res_flow = 0;
        for (int &i = iter[v]; i < (int)G[v].size(); ++i) {
            FlowEdge<FLOW> &e = G[v][i], &re = G.get_rev_edge(e);
            if (level[v] >= level[e.to] || e.cap == 0) continue;
            FLOW flow = self(self, e.to, min(up_flow - res_flow, e.cap));
            if (flow <= 0) continue;
            res_flow += flow;
            e.cap -= flow, e.flow += flow;
            re.cap += flow, re.flow -= flow;
            if (res_flow == up_flow) break;
        }
        return res_flow;
    };
    
    // flow
    while (current_flow < limit_flow) {
        bfs();
        if (level[t] < 0) break;
        iter.assign((int)iter.size(), 0);
        while (current_flow < limit_flow) {
            FLOW flow = dfs(dfs, s, limit_flow - current_flow);
            if (!flow) break;
            current_flow += flow;
        }
    }
    return current_flow;
};

template<class FLOW> FLOW Dinic(FlowGraph<FLOW> &G, int s, int t) {
    return Dinic(G, s, t, numeric_limits<FLOW>::max());
}


// edge class (for min-cost flow)
template<class FLOW, class COST> struct FlowCostEdge {
    // core members
    int rev, from, to;
    FLOW cap, icap, flow;
    COST cost;
    
    // constructor
    FlowCostEdge() {}
    FlowCostEdge(int rev, int from, int to, FLOW cap, COST cost)
    : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(0), cost(cost) {}
    void reset() { cap = icap, flow = 0; }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowCostEdge& e) {
        return s << e.from << " -> " << e.to << " (" << e.flow << "/" << e.icap << ", " << e.cost << ")";
    }
};

// graph class (for min-cost flow)
template<class FLOW, class COST> struct FlowCostGraph {
    // core members
    vector<vector<FlowCostEdge<FLOW, COST>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    vector<COST> pot; // pot[v] := potential (e.cost + pot[e.from] - pos[e.to] >= 0)
    bool include_negative_edge = false;
    
    // constructor
    FlowCostGraph(int n = 0) : list(n), pot(n), include_negative_edge(false) { }
    void init(int n = 0) {
        list.clear(), list.resize(n);
        pos.clear();
        pot.assign(n, 0);
        include_negative_edge = false;
    }
    
    // getter
    vector<FlowCostEdge<FLOW, COST>> &operator [] (int i) {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    const vector<FlowCostEdge<FLOW, COST>> &operator [] (int i) const {
        assert(0 <= i && i < list.size());
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowCostEdge<FLOW, COST> &get_rev_edge(const FlowCostEdge<FLOW, COST> &e) {
        return list[e.to][e.rev];
    }
    FlowCostEdge<FLOW, COST> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowCostEdge<FLOW, COST> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowCostEdge<FLOW, COST>> get_edges() const {
        vector<FlowCostEdge<FLOW, COST>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowCostEdge<FLOW, COST> &e : list[i]) e.reset();
        }
    }
    
    // add_edge
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, 0, -cost));
        if (cost < 0) include_negative_edge = true;
    }
    void add_edge(int from, int to, FLOW cap, FLOW rcap, COST cost) {
        assert(0 <= from && from < list.size() && 0 <= to && to < list.size());
        assert(cap >= 0);
        int from_id = int(list[from].size()), to_id = int(list[to].size());
        if (from == to) to_id++;
        pos.emplace_back(from, from_id);
        list[from].push_back(FlowCostEdge<FLOW, COST>(to_id, from, to, cap, cost));
        list[to].push_back(FlowCostEdge<FLOW, COST>(from_id, to, from, rcap, -cost));
        if (cost < 0) include_negative_edge = true;
    }

    // find initial potential (to resolve initial negative-edge)
    // pot[v] := potential (e.cost + pot[e.from] - pos[e.to] >= 0)
    bool calc_potential_dag() {
        pot.assign(size(), 0);
        vector<int> deg(size(), 0), st;
        for (int v = 0; v < size(); v++) for (const auto &e : list[v]) deg[e.to] += (e.cap > 0);
        st.reserve(size());
        for (int v = 0; v < size(); v++) if (!deg[v]) st.emplace_back(v);
        for (int i = 0; i < size(); i++) {
            if (st.size() == i) return false;  // not DAG
            int cur = st[i];
            for (const auto &e : list[cur]) {
                if (!e.cap) continue;
                deg[e.to]--;
                if (deg[e.to] == 0) st.emplace_back(e.to);
                if (pot[e.to] >= pot[cur] + e.cost) pot[e.to] = pot[cur] + e.cost;
            }
        }
        return true;
    }
    bool calc_potential_spfa() {
        pot.assign(size(), 0);
        queue<int> que;
        vector<bool> inque(size(), false);
        vector<int> cnt(size(), 0);
        for (int v = 0; v < size(); v++) que.push(v), inque[v] = true;
        while (!que.empty()) {
            int cur = que.front();
            que.pop();
            inque[cur] = false;
            if (cnt[cur] > size()) return false;  // include negative-cycle
            cnt[cur]++;
            for (const auto &e : list[cur]) {
                if (!e.cap) continue;
                if (pot[e.to] > pot[cur] + e.cost) {
                    pot[e.to] = pot[cur] + e.cost;
                    if (!inque[e.to]) inque[e.to] = true, que.push(e.to);
                }
            }
        }
        return true;
    }
    bool calc_potential() {
        return calc_potential_dag() || calc_potential_spfa();
    }
    bool init_potential() {
        if (!include_negative_edge) return true;
        return calc_potential();
    }

    // debug
    friend ostream& operator << (ostream& s, const FlowCostGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// min-cost max-flow (<= limit_flow), slope ver.
template<class FLOW, class COST> vector<pair<FLOW, COST>>
MinCostFlowSlope(FlowCostGraph<FLOW, COST> &G, int S, int T, FLOW limit_flow)
{
    // result values
    FLOW cur_flow = 0;
    COST cur_cost = 0, pre_cost = numeric_limits<COST>::max() / 2;
    vector<pair<FLOW, COST>> res;
    res.emplace_back(cur_flow, cur_cost);
    
    // intermediate values
    vector<COST> dist((int)G.size(), numeric_limits<COST>::max() / 2);
    vector<int> prevv((int)G.size(), -1), preve((int)G.size(), -1);
    
    // dual
    auto dual_step = [&]() -> bool {
        dist.assign((int)G.size(), numeric_limits<COST>::max() / 2);
        dist[S] = 0;
        priority_queue<pair<COST,int>, vector<pair<COST,int>>, greater<pair<COST,int>>> que;
        que.emplace(0, S);
        while (!que.empty()) {
            auto [cur, v] = que.top();
            que.pop();
            if (cur > dist[v]) continue;
            for (int i = 0; i < (int)G[v].size(); i++) {
                const auto &e = G[v][i];
                COST add = e.cost + G.pot[v] - G.pot[e.to];
                if (e.cap && dist[e.to] > dist[v] + add) {
                    dist[e.to] = dist[v] + add;
                    prevv[e.to] = v;
                    preve[e.to] = i;
                    que.emplace(dist[e.to], e.to);
                }
            }
        }
        return dist[T] < numeric_limits<COST>::max() / 2;
    };
    
    // primal
    auto primal_step = [&]() -> void {
        for (int v = 0; v < G.size(); v++) {
            if (dist[v] < numeric_limits<COST>::max() / 2) G.pot[v] += dist[v];
            else G.pot[v] = numeric_limits<COST>::max() / 2;
        }
        FLOW flow = limit_flow - cur_flow;
        COST cost = G.pot[T] - G.pot[S];
        for (int v = T; v != S; v = prevv[v]) {
            flow = min(flow, G[prevv[v]][preve[v]].cap);
        }
        for (int v = T; v != S; v = prevv[v]) {
            FlowCostEdge<FLOW, COST> &e = G[prevv[v]][preve[v]];
            FlowCostEdge<FLOW, COST> &re = G.get_rev_edge(e);
            e.cap -= flow, e.flow += flow;
            re.cap += flow, re.flow -= flow;
        }
        cur_flow += flow;
        cur_cost += flow * cost;
        if (pre_cost == cost) res.pop_back();
        res.emplace_back(cur_flow, cur_cost);
        pre_cost = cost;
    };

    // initialize potential
    assert(G.init_potential());
    
    // primal-dual
    while (cur_flow < limit_flow) {
        if (!dual_step()) break;
        primal_step();
    }
    return res;
}

// min-cost max-flow, slope ver.
template<class FLOW, class COST> vector<pair<FLOW, COST>>
MinCostFlowSlope(FlowCostGraph<FLOW, COST> &G, int S, int T)
{
    return MinCostFlowSlope(G, S, T, numeric_limits<FLOW>::max());
}

// min-cost max-flow (<= limit_flow)
template<class FLOW, class COST> pair<FLOW, COST>
MinCostFlow(FlowCostGraph<FLOW, COST> &G, int S, int T, FLOW limit_flow)
{
    return MinCostFlowSlope(G, S, T, limit_flow).back();
}

// min-cost max-flow (<= limit_flow)
template<class FLOW, class COST> pair<FLOW, COST>
MinCostFlow(FlowCostGraph<FLOW, COST> &G, int S, int T)
{
    return MinCostFlow(G, S, T, numeric_limits<FLOW>::max());
}


// Min Cost Circulation Flow by Cost-Scaling 
template<class FLOW, class COST> COST MinCostCirculation(FlowCostGraph<FLOW, COST> &G) {
    COST eps = 0;
    vector<FLOW> balance(G.size(), 0);
    vector<COST> price(G.size(), 0);
    
    auto newcost = [&](const FlowCostEdge<FLOW, COST> &e) -> COST {
        return e.cost * (COST)G.size() - price[e.from] + price[e.to];
    };

    auto ConstructGaux = [&]() -> void {
        vector<bool> visited(G.size(), false);
        auto dfs = [&](auto &&dfs, int v) -> void {
            visited[v] = true;
            for (int i = 0; i < G[v].size(); ++i) {
                FlowCostEdge<FLOW, COST> &e = G[v][i];
                if (e.cap > 0 && !visited[e.to] && newcost(e) < 0) dfs(dfs, e.to);
            }
        };
        for (int v = 0; v < G.size(); ++v) if (balance[v] > 0) dfs(dfs, v);
        for (int v = 0; v < G.size(); ++v) if (visited[v]) price[v] += eps;
    };

    auto augment_blocking_flow = [&]() -> bool {
        vector<COST> iter(G.size(), 0);
        auto augment = [&](auto &&augment, int v, FLOW flow) -> FLOW {
            if (balance[v] < 0) {
                FLOW dif = min(flow, -balance[v]);
                balance[v] += dif;
                return dif;
            }
            for (; iter[v] < G[v].size(); iter[v]++) {
                auto &e = G[v][iter[v]], &re = G.get_rev_edge(e);
                if (e.cap > 0 && newcost(e) < 0) {
                    FLOW dif = augment(augment, e.to, min(flow, e.cap));
                    if (dif > 0) {
                        e.cap -= dif, e.flow += dif;
                        re.cap += dif, re.flow -= dif;
                        return dif;
                    }
                }
            }
            return 0;
        };
        bool finish = true;
        for (int v = 0; v < G.size(); ++v) {
            FLOW flow;
            while (balance[v] > 0 && (flow = augment(augment, v, balance[v])) > 0)
                balance[v] -= flow;
            if (balance[v] > 0) finish = false;
        }
        if (finish) return true;
        else return false;
    };

    for (int v = 0; v < G.size(); ++v) {
        for (int i = 0; i < G[v].size(); ++i) {
            FlowCostEdge<FLOW, COST> &e = G[v][i];
            if (e.cap > 0) eps = max(eps, -e.cost * (COST)G.size());
        }
        price[v] = 0;
    }
    while (eps > 1) {
        eps /= 2;
        for (int v = 0; v < G.size(); ++v) {
            for (int i = 0; i < G[v].size(); ++i) {
                auto &e = G[v][i], &re = G.get_rev_edge(e);
                if (e.cap > 0 && newcost(e) < 0) {
                    FLOW flow = e.cap;
                    balance[e.from] -= flow, balance[e.to] += flow;
                    e.cap -= flow, e.flow += flow;
                    re.cap += flow, re.flow -= flow;
                }
            }
        }
        while (true) {
            ConstructGaux();
            if (augment_blocking_flow()) break;
        }
    }
    COST res = 0;
    const auto &edges = G.get_edges();
    for (const auto &e : edges) res += e.flow * e.cost;
    return res;
}

// Minimum Cost b-flow (come down to min-cost circulation)
template<class FLOW, class COST> struct MinCostBFlow {
    // Edge
    struct Edge {
        int from, to;
        FLOW lower_cap, upper_cap, flow;
        COST cost;
    };

    // inner values
    int V;
    vector<Edge> edges;
    vector<FLOW> dss;  // demand (< 0) and supply (> 0)
    vector<COST> dual;
    FlowCostGraph<FLOW, COST> G;

    // constructor
    MinCostBFlow() {}
    MinCostBFlow(int V) : V(V), dss(V) {}

    // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        edges.push_back({from, to, 0, cap, 0, cost});
    }
    void add_edge(int from, int to, FLOW lower_cap, FLOW upper_cap, COST cost) {
        assert(lower_cap <= upper_cap);
        edges.push_back({from, to, lower_cap, upper_cap, 0, cost});
    }
    void add_ds(int v, FLOW ds) {
        assert(0 <= v && v < V);
        dss[v] += ds;
    }

    // solver
    pair<bool, COST> solve(bool calc_potential = true) {
        // feasibility check
        FlowGraph<FLOW> sg(V + 2);
        int s = V, t = s + 1;
        for (const auto &e : edges) {
            dss[e.to] += e.lower_cap, dss[e.from] -= e.lower_cap;
            sg.add_edge(e.from, e.to, e.upper_cap - e.lower_cap);
        }
        FLOW ssum = 0, tsum = 0;
        for (int i = 0; i < V; i++) {
            if (dss[i] > 0) ssum += dss[i], sg.add_edge(s, i, dss[i]);
            else if (dss[i] < 0) tsum -= dss[i], sg.add_edge(i, t, -dss[i]);
        }
        if (ssum != tsum) return {false, COST(0)};
        if (Dinic(sg, s, t) < ssum) return {false, COST(0)};

        // come down to min-cost circulation
        G.init(V);
        for (int i = 0; i < (int)edges.size(); i++) {
            auto &e = edges[i];
            const auto &ge = sg.get_edge(i);
            G.add_edge(ge.from, ge.to, ge.cap, ge.flow, e.cost);
        }
        MinCostCirculation(G);

        // find min-cost
        COST res = 0;
        for (int i = 0; i < (int)edges.size(); i++) {
            auto &e = edges[i];
            const auto &ge = G.get_edge(i);
            e.flow = e.upper_cap - ge.cap;
            res += e.flow * e.cost;
        }

        // find dual
        if (calc_potential) {
            G.calc_potential();
            dual = G.pot;
        }
        return {true, res};
    }
};

// upper and lower demand-supply
template<class FLOW, class COST> struct GeneraliezedMinCostBFlow {
    // Edge
    struct Edge {
        int from, to;
        FLOW lower_cap, upper_cap, flow;
        COST cost;
    };

    // inner values
    int V;
    vector<Edge> edges;
    vector<FLOW> lower_dss, upper_dss;
    MinCostBFlow<FLOW, COST> mcbf;

    // constructor
    GeneraliezedMinCostBFlow() {}
    GeneraliezedMinCostBFlow(int V) : V(V), lower_dss(V), upper_dss(V), mcbf(V) {}

    // setter
        // setter
    void add_edge(int from, int to, FLOW cap, COST cost) {
        assert(cap >= 0);
        edges.push_back({from, to, 0, cap, 0, cost});
    }
    void add_edge(int from, int to, FLOW lower_cap, FLOW upper_cap, COST cost) {
        assert(lower_cap <= upper_cap);
        edges.push_back({from, to, lower_cap, upper_cap, 0, cost});
    }
    void add_ds(int v, FLOW lower_ds, FLOW upper_ds) {
        assert(0 <= v && v < V);
        lower_dss[v] = lower_ds, upper_dss[v] = upper_ds;
    }

    // solver
    pair<bool, COST> solve() {
        
    }
};

// Hopcroft-Karp
struct HopcroftKarp {
    const int NOT_MATCHED = -1;
    
    // input
    int size_left, size_right;
    vector<vector<int>> list; // left to right
    vector<vector<int>> rlist; // right to left

    // results
    vector<int> lr, rl;
    
    // intermediate results
    vector<bool> seen, matched;
    vector<int> level;
    
    // constructor
    HopcroftKarp(int L, int R) : size_left(L), size_right(R), list(L), rlist(R) {}
    void add_edge(int from, int to) {
        assert(from >= 0 && from < size_left);
        assert(to >= 0 && to < size_right);
        list[from].emplace_back(to);
        rlist[to].emplace_back(from);
    }

    // getter, debugger
    const vector<int> &operator [] (int i) const { 
        return list[i];
    }
    friend ostream& operator << (ostream& s, const HopcroftKarp& G) {
        s << endl;
        for (int i = 0; i < (int)G.list.size(); ++i) {
            s << i << ": ";
            for (int j = 0; j < (int)G.list[i].size(); ++j) {
                s << G.list[i][j];
                if (j + 1 != (int)G.list[i].size()) s << ", ";
            }
            s << endl;
        }
        return s;
    }
    
    // solver
    void hobfs() {
        queue<int> que;
        for (int left = 0; left < size_left; ++left) {
            level[left] = -1;
            if (!matched[left]) {
                que.push(left);
                level[left] = 0;
            }
        }
        level[size_left] = size_left;
        while (!que.empty()) {
            int left = que.front();
            que.pop();
            for (int right : list[left]) {
                int next = rl[right];
                if (next == NOT_MATCHED) next = size_left;
                if (level[next] == -1) {
                    level[next] = level[left] + 1;
                    que.push(next);
                }
            }
        }
    }
    bool hodfs(int left) {
        if (left == size_left) return true;
        if (seen[left]) return false;
        seen[left] = true;
        for (int right : list[left]) {
            int next = rl[right];
            if (next == NOT_MATCHED) next = size_left;
            if (level[next] > level[left] && hodfs(next)) {
                lr[left] = right;
                rl[right] = left;
                return true;
            }
        }
        return false;
    }
    int solve() {
        seen.assign(size_left, false);
        matched.assign(size_left, false);
        level.assign(size_left + 1, -1);
        lr.assign(size_left, -1);
        rl.assign(size_right, -1);
        int res = 0;
        while (true) {
            hobfs();
            seen.assign(size_left, false);
            bool finished = true;
            for (int left = 0; left < size_left; ++left) {
                if (!matched[left] && hodfs(left)) {
                    matched[left] = true;
                    ++res;
                    finished = false;
                }
            }
            if (finished) break;
        }
        for (int r = 0; r < size_right; r++) {
            if (rl[r] != NOT_MATCHED) lr[rl[r]] = r;
        }
        return res;
    }

    // various construction
    // max matching
    vector<pair<int,int>> get_matching() {
        vector<pair<int,int>> res;
        for (int v = 0; v < size_left; v++) {
            if (lr[v] == NOT_MATCHED) continue;
            res.emplace_back(v, lr[v]);
        }
        return res;
    }

    // enumerate reachable nodes (0: left, 1: right)
    const int LEFT = 0, RIGHT = 1;
    pair<vector<bool>, vector<bool>> get_reachable() {
        vector<bool> can_left(size_left, false);
        vector<bool> can_right(size_right, false);
        queue<pair<int,int>> que;
        for (int v = 0; v < size_left; v++) {
            if (lr[v] == NOT_MATCHED) {
                can_left[v] = true;
                que.push({LEFT, v});
            }
        }
        while (!que.empty()) {
            auto [which, v] = que.front();
            que.pop();
            if (which == LEFT) {
                for (auto r : list[v]) {
                    if (!can_right[r]) {
                        can_right[r] = true;
                        que.push({RIGHT, r});
                    }
                }
            } else {
                int l = rl[v];
                if (l != NOT_MATCHED && !can_left[l]) {
                    can_left[l] = true;
                    que.push({LEFT, l});
                }
            }
        }
        return {can_left, can_right};
    }

    // max independent set (0: left, 1: right)
    vector<pair<int,int>> get_independent_set() {
        vector<pair<int,int>> res;
        auto [can_left, can_right] = get_reachable();
        for (int v = 0; v < size_left; v++) {
            if (can_left[v]) res.emplace_back(LEFT, v);
        }
        for (int v = 0; v < size_right; v++) {
            if (!can_right[v]) res.emplace_back(RIGHT, v);
        }
        return res;
    }

    // min vertex-cover (0: left, 1: right)
    vector<pair<int,int>> get_vertex_cover() {
        vector<pair<int,int>> res;
        auto [can_left, can_right] = get_reachable();
        for (int v = 0; v < size_left; v++) {
            if (!can_left[v]) res.emplace_back(LEFT, v);
        }
        for (int v = 0; v < size_right; v++) {
            if (can_right[v]) res.emplace_back(RIGHT, v);
        }
        return res;
    }

    // min edge-cover (0: left, 1: right)
    vector<pair<int,int>> get_edge_cover() {
        vector<pair<int,int>> res = get_matching();
        for (int v = 0; v < size_left; v++) {
            if (list[v].empty()) return vector<pair<int,int>>();  // infeasible
            if (lr[v] == NOT_MATCHED) res.emplace_back(v, list[v][0]);
        }
        for (int v = 0; v < size_right; v++) {
            if (rlist[v].empty()) return vector<pair<int,int>>();  // infeasible
            if (rl[v] == NOT_MATCHED) res.emplace_back(rlist[v][0], v);
        }
        return res;
    }
};

// 1, 2, 3-variable submodular optimization
/*
 N 個の bool 変数 x_0, x_1, ..., x_{N-1} について、以下の形のコストが定められたときの最小コストを求める
 
 ・1 変数 xi に関するコスト (1 変数劣モジュラ関数)
    xi = F のときのコスト, xi = T のときのコスト
 
 ・2 変数 xi, xj 間の関係性についてのコスト (2 変数劣モジュラ関数)
 　　(xi, xj) = (F, F): コスト A
 　　(xi, xj) = (F, T): コスト B
 　　(xi, xj) = (T, F): コスト C
 　　(xi, xj) = (T, T): コスト D
 　(ただし、B + C >= A + D でなければならない)
 
 ・よくある例は、A = B = D = 0, C >= 0 の形である (特に関数化している)
    ・この場合は、特に Project Selection Problem と呼ばれ、俗に「燃やす埋める」などとも呼ばれる
    ・xi = T, xj = F のときにコスト C がかかる
 
 ・他に面白い例として、A = B = C = 0, D <= 0 の形もある (これも関数化している)
    ・xi = T, xj = T のときに (-D) の利得が得られる
 
 ・3 変数 xi, xj, xk 間の関係性についてのコスト (3 変数劣モジュラ関数)
 　　(xi, xj, xk) = (F, F, F): コスト A
 　　(xi, xj, xk) = (F, F, T): コスト B
 　　(xi, xj, xk) = (F, T, F): コスト C
 　　(xi, xj, xk) = (F, T, T): コスト D
 　　(xi, xj, xk) = (T, F, F): コスト E
 　　(xi, xj, xk) = (T, F, T): コスト F
 　　(xi, xj, xk) = (T, T, F): コスト G
 　　(xi, xj, xk) = (T, T, T): コスト H
 */
template<class COST> struct ThreeVariableSubmodularOpt {
    // constructors
    ThreeVariableSubmodularOpt() : N(2), S(0), T(0), OFFSET(0) {}
    ThreeVariableSubmodularOpt(int n, COST inf = numeric_limits<COST>::max() / 2)
    : N(n), S(n), T(n + 1), OFFSET(0), INF(inf), list(n + 2) {}
    
    // initializer
    void init(int n, COST inf = numeric_limits<COST>::max() / 2) {
        N = n, S = n, T = n + 1;
        OFFSET = 0, INF = inf;
        list.assign(N + 2, Edge());
        pos.clear();
    }

    // add 1-Variable submodular functioin
    void add_single_cost(int xi, COST false_cost, COST true_cost) {
        assert(0 <= xi && xi < N);
        if (false_cost >= true_cost) {
            OFFSET += true_cost;
            add_edge(S, xi, false_cost - true_cost);
        } else {
            OFFSET += false_cost;
            add_edge(xi, T, true_cost - false_cost);
        }
    }
    
    // add "project selection" constraint
    // xi = T, xj = F: strictly prohibited
    void add_psp_constraint(int xi, int xj) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        add_edge(xi, xj, INF);
    }
    
    // add "project selection" penalty
    // xi = T, xj = F: cost C
    void add_psp_penalty(int xi, int xj, COST C) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(C >= 0);
        add_edge(xi, xj, C);
    }
    
    // add both True profit
    // xi = T, xj = T: profit P (cost -P)
    void add_both_true_profit(int xi, int xj, COST P) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(P >= 0);
        OFFSET -= P;
        add_edge(S, xi, P);
        add_edge(xi, xj, P);
    }
    
    // add both False profit
    // xi = F, xj = F: profit P (cost -P)
    void add_both_false_profit(int xi, int xj, COST P) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(P >= 0);
        OFFSET -= P;
        add_edge(xj, T, P);
        add_edge(xi, xj, P);
    }
    
    // add general 2-variable submodular function
    // (xi, xj) = (F, F): A, (F, T): B
    // (xi, xj) = (T, F): C, (T, T): D
    void add_submodular_function(int xi, int xj, COST A, COST B, COST C, COST D) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(B + C >= A + D);  // assure submodular function
        OFFSET += A;
        add_single_cost(xi, 0, D - B);
        add_single_cost(xj, 0, B - A);
        add_psp_penalty(xi, xj, B + C - A - D);
    }
    
    // add all True profit
    // y = F: not gain profit (= cost is P), T: gain profit (= cost is 0)
    // y: T, xi: F is prohibited
    void add_all_true_profit(const vector<int> &xs, COST P) {
        assert(P >= 0);
        int y = (int)list.size();
        list.resize(y + 1);
        OFFSET -= P;
        add_edge(S, y, P);
        for (auto xi : xs) {
            assert(xi >= 0 && xi < N);
            add_edge(y, xi, INF);
        }
    }
    
    // add all False profit
    // y = F: gain profit (= cost is 0), T: not gain profit (= cost is P)
    // xi = T, y = F is prohibited
    void add_all_false_profit(const vector<int> &xs, COST P) {
        assert(P >= 0);
        int y = (int)list.size();
        list.resize(y + 1);
        OFFSET -= P;
        add_edge(y, T, P);
        for (auto xi : xs) {
            assert(xi >= 0 && xi < N);
            add_edge(xi, y, INF);
        }
    }
    
    // add general 3-variable submodular function
    // (xi, xj, xk) = (F, F, F): cost A
    // (xi, xj, xk) = (F, F, T): cost B
    // (xi, xj, xk) = (F, T, F): cost C
    // (xi, xj, xk) = (F, T, T): cost D
    // (xi, xj, xk) = (T, F, F): cost E
    // (xi, xj, xk) = (T, F, T): cost F
    // (xi, xj, xk) = (T, T, F): cost G
    // (xi, xj, xk) = (T, T, T): cost H
    void add_submodular_function(int xi, int xj, int xk,
                                 COST A, COST B, COST C, COST D,
                                 COST E, COST F, COST G, COST H) {
        assert(0 <= xi && xi < N);
        assert(0 <= xj && xj < N);
        assert(0 <= xk && xk < N);
        COST P = (A + D + F + G) - (B + C + E + H);
        COST P12 = (C + E) - (A + G), P13 = (D + G) - (C + H);
        COST P21 = (D + F) - (B + H), P23 = (B + C) - (A + D);
        COST P31 = (B + E) - (A + F), P32 = (F + G) - (E + H);
        assert(P12 >= 0 && P21 >= 0);
        assert(P23 >= 0 && P32 >= 0);
        assert(P31 >= 0 && P13 >= 0);
        if (P >= 0) {
            OFFSET += A;
            add_single_cost(xi, 0, F - B);
            add_single_cost(xj, 0, G - E);
            add_single_cost(xk, 0, D - C);
            add_psp_penalty(xj, xi, P12);
            add_psp_penalty(xk, xj, P23);
            add_psp_penalty(xi, xk, P31);
            add_all_true_profit({xi, xj, xk}, P);
        } else {
            OFFSET += H;
            add_single_cost(xi, C - G, 0);
            add_single_cost(xj, B - D, 0);
            add_single_cost(xk, E - F, 0);
            add_psp_penalty(xi, xj, P21);
            add_psp_penalty(xj, xk, P32);
            add_psp_penalty(xk, xi, P13);
            add_all_false_profit({xi, xj, xk}, -P);
        }
    }
    
    // solve
    COST solve() {
        return dinic() + OFFSET;
    }
    
    // reconstrcut the optimal assignment
    vector<bool> reconstruct() {
        vector<bool> res(N, false), seen(list.size(), false);
        queue<int> que;
        seen[S] = true;
        que.push(S);
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            for (const auto &e : list[v]) {
                if (e.cap && !seen[e.to]) {
                    if (e.to < N) res[e.to] = true;
                    seen[e.to] = true;
                    que.push(e.to);
                }
            }
        }
        return res;
    }
    
    // debug
    friend ostream& operator << (ostream& s, const ThreeVariableSubmodularOpt &tvs) {
        const auto &edges = tvs.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
    
private:
    // edge class
    struct Edge {
        // core members
        int rev, from, to;
        COST cap, icap, flow;
        
        // constructor
        Edge(int r, int f, int t, COST c)
        : rev(r), from(f), to(t), cap(c), icap(c), flow(0) {}
        void reset() { cap = icap, flow = 0; }
        
        // debug
        friend ostream& operator << (ostream& s, const Edge& E) {
            return s << E.from << "->" << E.to << '(' << E.flow << '/' << E.icap << ')';
        }
    };
    
    // inner data
    int N, S, T;
    COST OFFSET, INF;
    vector<vector<Edge>> list;
    vector<pair<int,int>> pos;
    
    // add edge
    Edge &get_rev_edge(const Edge &e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    Edge &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const Edge &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<Edge> get_edges() const {
        vector<Edge> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    void add_edge(int from, int to, COST cap) {
        if (!cap) return;
        pos.emplace_back(from, (int)list[from].size());
        list[from].push_back(Edge((int)list[to].size(), from, to, cap));
        list[to].push_back(Edge((int)list[from].size() - 1, to, from, 0));
    }
    
    // Dinic's algorithm
    COST dinic(COST limit_flow) {
        COST current_flow = 0;
        vector<int> level((int)list.size(), -1), iter((int)list.size(), 0);
        
        // Dinic BFS
        auto bfs = [&]() -> void {
            level.assign((int)list.size(), -1);
            level[S] = 0;
            queue<int> que;
            que.push(S);
            while (!que.empty()) {
                int v = que.front();
                que.pop();
                for (const Edge &e : list[v]) {
                    if (level[e.to] < 0 && e.cap > 0) {
                        level[e.to] = level[v] + 1;
                        if (e.to == T) return;
                        que.push(e.to);
                    }
                }
            }
        };
        
        // Dinic DFS
        auto dfs = [&](auto self, int v, COST up_flow) {
            if (v == T) return up_flow;
            COST res_flow = 0;
            for (int &i = iter[v]; i < (int)list[v].size(); ++i) {
                Edge &e = list[v][i], &re = get_rev_edge(e);
                if (level[v] >= level[e.to] || e.cap == 0) continue;
                COST flow = self(self, e.to, min(up_flow - res_flow, e.cap));
                if (flow <= 0) continue;
                res_flow += flow;
                e.cap -= flow, e.flow += flow;
                re.cap += flow, re.flow -= flow;
                if (res_flow == up_flow) break;
            }
            return res_flow;
        };
        
        // flow
        while (current_flow < limit_flow) {
            bfs();
            if (level[T] < 0) break;
            iter.assign((int)iter.size(), 0);
            while (current_flow < limit_flow) {
                COST flow = dfs(dfs, S, limit_flow - current_flow);
                if (!flow) break;
                current_flow += flow;
            }
        }
        return current_flow;
    };
    COST dinic() {
        return dinic(numeric_limits<COST>::max() / 2);
    }
};


//------------------------------//
// Union-Find
//------------------------------//

// Union-Find
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
    
    bool merge(int x, int y, bool merge_technique = false) {
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


//------------------------------//
// Sparse Table, Binary Indexed Tree
//------------------------------//

// BIT
template <class Abel> struct BIT {
    Abel UNITY_SUM = 0;
    vector<Abel> dat;
    
    // [0, n)
    BIT(int n, Abel unity = 0) : UNITY_SUM(unity), dat(n, unity) { }
    void init(int n) {
        dat.assign(n, UNITY_SUM);
    }
    int size() const {
        return (int)dat.size();
    }
    
    // a is 0-indexed
    inline void add(int a, Abel x) {
        for (int i = a; i < (int)dat.size(); i |= i + 1)
            dat[i] = dat[i] + x;
    }
    
    // [0, a), a is 0-indexed, [a, b), a and b are 0-indexed
    inline Abel sum(int a) const {
        Abel res = UNITY_SUM;
        for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[i];
        return res;
    }
    inline Abel sum(int a, int b) const {
        return sum(b) - sum(a);
    }
    inline Abel operator [] (int i) const {
        return sum(i, i + 1);
    }
    
    // debug
    friend ostream& operator << (ostream &s, const BIT &bit) {
        for (int i = 0; i < (int)bit.size(); ++i) s << bit[i] << " ";
        return s;
    }
};

// RangeAddRangeSum by BIT
template<class Abel> struct RangeAddRangeSum {
    Abel UNITY_SUM = 0;
    vector<Abel> dat[2];

    // [0, n)
    RangeAddRangeSum(int n, Abel unity = 0) : UNITY_SUM(unity) {
        init(n);
    }
    void init(int n) {
        for (int iter = 0; iter < 2; ++iter)
            dat[iter].assign(n + 1, UNITY_SUM);
    }
    int size() const {
        return (int)dat[0].size();
    }
    
    // [a, b), a and b are 0-indexed
    inline void sub_add(int p, int a, Abel x) {
        for (int i = a; i < (int)dat[p].size(); i |= i + 1)
            dat[p][i] = dat[p][i] + x;
    }
    inline void add(int a, int b, Abel x) {
        sub_add(0, a, x * (-a));
        sub_add(1, a, x);
        sub_add(0, b, x * b);
        sub_add(1, b, x * (-1));
    }
    
    // [a, b), a and b are 0-indexed
    inline Abel sub_sum(int p, int a) {
        Abel res = UNITY_SUM;
        for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[p][i];
        return res;
    }
    inline Abel sum(int a, int b) {
        return sub_sum(0, b)
            + sub_sum(1, b) * b
            - sub_sum(0, a)
            - sub_sum(1, a) * a;
    }
    inline Abel operator [] (int i) const {
        return sum(i, i + 1);
    }
    
    // debug
    friend ostream& operator << (ostream &s, const RangeAddRangeSum &bit) {
        for (int i = 0; i < (int)bit.size(); ++i) s << bit[i] << " ";
        return s;
    }
};

// Fast MultiSet By BIT
template<class Abel = int> struct FastMultiSetByBIT {
    int topbit(int x) const { return (x == 0 ? -1 : 31 - __builtin_clz(x)); }
    int lowbit(int x) const { return (x == 0 ? -1 : __builtin_ctz(x)); }
    int lim;
    Abel IDENTITY;
    vector<Abel> dat;
    
    // [0, n)
    FastMultiSetByBIT(int n, Abel identity = 0)
    : lim(n), IDENTITY(identity), dat(n, identity) { }
    void init(int n, Abel identity = 0) {
        lim = n;
        IDENTITY = identity;
        dat.assign(n, IDENTITY);
    }
    
    // p is 0-indexed
    void add(int p, Abel x) {
        if (p < 0) p = 0;
        for (int i = p; i < (int)dat.size(); i |= i + 1)
            dat[i] = dat[i] + x;
    }
    
    // [0, p), p is 0-indexed
    Abel sum(int p) const {
        if (p > lim) p = lim;
        Abel res = IDENTITY;
        for (int i = p - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[i];
        return res;
    }
    
    // [l, r), l and r are 0-indexed
    Abel sum(int l, int r) const {
        return sum(r) - sum(l);
    }
    
    // insert, erase, count, min, max
    void insert(int x, Abel num = 1) { add(x, num); }
    void erase(int x, Abel num = 1) { add(x, -min(num, count(x))); }
    void clear() { dat.assign(lim, IDENTITY); }
    Abel count(int x) const { return sum(x, x + 1); }
    Abel count(int l, int r) const { return sum(l, r); }
    Abel size() const { return sum(lim); }
    int get_min() const { return next(); }
    int get_max() const { return prev(); }
    
    // get max r s.t. check(sum(l, r)) = True (0-indexed), O(log N)
    // check(IDENTITY) must be True
    int max_right(const function<bool(Abel)> check, int l = 0) const {
        if (l >= lim) return lim;
        assert(check(IDENTITY));
        Abel s = IDENTITY;
        int k = 0;
        while (true) {
            if (l % 2 == 1) s = s - dat[l - 1], --l;
            if (l <= 0) {
                k = topbit(lim) + 1;
                break;
            }
            k = lowbit(l) - 1;
            if (l + (1 << k) > lim) break;
            if (!check(s + dat[l + (1 << k) - 1])) break;
            s = s - dat[l - 1];
            l -= l & -l;
        }
        while (k) {
            --k;
            if (l + (1 << k) - 1 < lim) {
                Abel ns = s + dat[l + (1 << k) - 1];
                if (check(ns)) {
                    l += (1 << k);
                    s = ns;
                }
            }
        }
        return l;
    }
    
    // get min l s.t. check(sum(l, r))  = True (0-indexed), O(log N)
    // check(IDENTITY) must be True
    int min_left(const function<bool(Abel)> check, int r = -1) const {
        if (r == -1) r = lim;
        if (r <= 0) return 0;
        assert(check(IDENTITY));
        Abel s = IDENTITY;
        int k = 0;
        while (r > 0 && check(s)) {
            s = s + dat[r - 1];
            k = lowbit(r);
            r -= r & -r;
        }
        if (check(s)) return 0;
        while (k) {
            --k;
            Abel ns = s - dat[r + (1 << k) - 1];
            if (!check(ns)) {
                r += (1 << k);
                s = ns;
            }
        }
        return r + 1;
    }
              
    // k-th number that is not less than l (k is 0-indexed)
    int get(Abel k, int l = 0) const {
        return max_right([&](Abel x) { return x <= k; }, l);
    }
    int operator [] (int k) const { return get(k); }
    
    // next (including x)
    int next(int l = 0) const {
        if (l < 0) l = 0;
        if (l > lim) l = lim;
        return max_right([&](Abel x) { return x <= 0; }, l);
    }
    
    // prev (including x)
    int prev(int r) const {
        if (r > lim) r = lim;
        return min_left([&](Abel x) { return x <= 0; }, r + 1) - 1;
    }
    int prev() const {
        return prev(lim);
    }
    
    // debug
    friend ostream& operator << (ostream &s, const FastMultiSetByBIT &fs) {
        for (int x = fs.get_min(); x < fs.lim; x = fs.next(x + 1)) {
            for (int i = 0; i < fs.count(x); i++) s << x << " ";
        }
        return s;
    }
};

// Sparse Table
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


//------------------------------//
// Segment Tree
//------------------------------//

// Segment Tree
template<class Monoid> struct SegmentTree {
    using Func = function<Monoid(Monoid, Monoid)>;

    // core member
    int N;
    Func OP;
    Monoid IDENTITY;
    
    // inner data
    int log, offset;
    vector<Monoid> dat;

    // constructor
    SegmentTree() {}
    SegmentTree(const Func &op, const Monoid &identity) : OP(op), IDENTITY(identity) { }
    SegmentTree(int n, const Func &op, const Monoid &identity) {
        init(n, op, identity);
    }
    SegmentTree(const vector<Monoid> &v, const Func &op, const Monoid &identity) {
        init(v, op, identity);
    }
    void init(const Func &op, const Monoid &identity) {
        OP = op;
        IDENTITY = identity;
    }
    void init(int n) {
        N = n;
        log = 0, offset = 1;
        while (offset < N) ++log, offset <<= 1;
        dat.assign(offset * 2, IDENTITY);
    }
    void init(const vector<Monoid> &v) {
        init((int)v.size());
        build(v);
    } 
    void init(int n, const Func &op, const Monoid &identity) {
        init(op, identity);
        init(n);
    }
    void init(const vector<Monoid> &v, const Func &op, const Monoid &identity) {
        init((int)v.size(), op, identity);
        build(v);
    }
    int size() const {
        return N;
    }

    // pull, build
    void pull(int k) {
        dat[k] = OP(dat[k * 2], dat[k * 2 + 1]);
    }
    void build(const vector<Monoid> &v) {
        assert(N == (int)v.size());
        for (int i = 0; i < N; ++i) dat[i + offset] = v[i];
        for (int k = offset - 1; k > 0; --k) pull(k);
    }
    
    // setter and getter, set: update A[i], i is 0-indexed, O(log N)
    void set(int i, const Monoid &v) {
        assert(0 <= i && i < N);
        int k = i + offset;
        dat[k] = v;
        while (k >>= 1) pull(k);
    }
    Monoid get(int i) const {
        assert(0 <= i && i < N);
        return dat[i + offset];
    }
    Monoid operator [] (int i) const {
        return get(i);
    }
    
    // get [l, r), l and r are 0-indexed, O(log N)
    Monoid prod(int l, int r) {
        assert(0 <= l && l <= r && r <= N);
        Monoid val_left = IDENTITY, val_right = IDENTITY;
        l += offset, r += offset;
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) val_left = OP(val_left, dat[l++]);
            if (r & 1) val_right = OP(dat[--r], val_right);
        }
        return OP(val_left, val_right);
    }
    Monoid all_prod() {
        return dat[1];
    }
    
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
    friend ostream& operator << (ostream &s, const SegmentTree &seg) {
        for (int i = 0; i < (int)seg.size(); ++i) {
            s << seg[i];
            if (i != (int)seg.size() - 1) s << " ";
        }
        return s;
    }

    // dump
    void dump() {
        int pt = 1;
        for (int h = 0; h <= log; h++) {
            for (int i = 0; i < (1<<h); i++) cout << dat[pt++] << " ";
            cout << endl;
        }
    }
};

// Lazy Segment Tree
template<class Monoid, class Action> struct LazySegmentTree {
    // various function types
    using FuncMonoid = function<Monoid(Monoid, Monoid)>;
    using FuncAction = function<Monoid(Action, Monoid)>;
    using FuncComposition = function<Action(Action, Action)>;

    // core member
    int N;
    FuncMonoid OP;
    FuncAction ACT;
    FuncComposition COMP;
    Monoid IDENTITY_MONOID;
    Action IDENTITY_ACTION;
    
    // inner data
    int log, offset;
    vector<Monoid> dat;
    vector<Action> lazy;
    
    // constructor
    LazySegmentTree() {}
    LazySegmentTree(const FuncMonoid op, const FuncAction act, const FuncComposition comp,
                    const Monoid &identity_monoid, const Action &identity_action) 
                    : OP(op), ACT(act), COMP(comp), 
                    IDENTITY_MONOID(identity_monoid), IDENTITY_ACTION(identity_action) {}
    LazySegmentTree(int n, const FuncMonoid op, const FuncAction act, const FuncComposition comp,
                    const Monoid &identity_monoid, const Action &identity_action) {
        init(n, op, act, comp, identity_monoid, identity_action);
    }
    LazySegmentTree(const vector<Monoid> &v,
                    const FuncMonoid op, const FuncAction act, const FuncComposition comp,
                    const Monoid &identity_monoid, const Action &identity_action) {
        init(v, op, act, comp, identity_monoid, identity_action);
    }
    void init(const FuncMonoid op, const FuncAction act, const FuncComposition comp,
              const Monoid &identity_monoid, const Action &identity_action) {
        OP = op, ACT = act, COMP = comp;
        IDENTITY_MONOID = identity_monoid, IDENTITY_ACTION = identity_action;      
    }
    void init(int n) {
        N = n, 
        log = 0, offset = 1;
        while (offset < N) ++log, offset <<= 1;
        dat.assign(offset * 2, IDENTITY_MONOID);
        lazy.assign(offset * 2, IDENTITY_ACTION);
    }
    void init(const vector<Monoid> &v) {
        init((int)v.size());
        build(v);
    }
    void init(int n, const FuncMonoid op, const FuncAction act, const FuncComposition comp,
              const Monoid &identity_monoid, const Action &identity_action) {
        init(op, act, comp, identity_monoid, identity_action);
        init(n);
    }
    void init(const vector<Monoid> &v,
              const FuncMonoid op, const FuncAction act, const FuncComposition comp,
              const Monoid &identity_monoid, const Action &identity_action) {
        init((int)v.size(), op, act, comp, identity_monoid, identity_action);
        build(v);
    }
    void build(const vector<Monoid> &v) {
        assert(N == (int)v.size());
        for (int i = 0; i < N; ++i) dat[i + offset] = v[i];
        for (int k = offset - 1; k > 0; --k) pull_dat(k);
    }
    int size() const {
        return N;
    }
    
    // basic functions for lazy segment tree
    void pull_dat(int k) {
        dat[k] = OP(dat[k * 2], dat[k * 2 + 1]);
    }
    void apply_lazy(int k, const Action &f) {
        dat[k] = ACT(f, dat[k]);
        if (k < offset) lazy[k] = COMP(f, lazy[k]);
    }
    void push_lazy(int k) {
        apply_lazy(k * 2, lazy[k]);
        apply_lazy(k * 2 + 1, lazy[k]);
        lazy[k] = IDENTITY_ACTION;
    }
    void pull_dat_deep(int k) {
        for (int h = 1; h <= log; ++h) pull_dat(k >> h);
    }
    void push_lazy_deep(int k) {
        for (int h = log; h >= 1; --h) push_lazy(k >> h);
    }
    
    // setter and getter, update A[i], i is 0-indexed, O(log N)
    void set(int i, const Monoid &v) {
        assert(0 <= i && i < N);
        int k = i + offset;
        push_lazy_deep(k);
        dat[k] = v;
        pull_dat_deep(k);
    }
    Monoid get(int i) {
        assert(0 <= i && i < N);
        int k = i + offset;
        push_lazy_deep(k);
        return dat[k];
    }
    Monoid operator [] (int i) {
        return get(i);
    }
    
    // apply f for index i
    void apply(int i, const Action &f) {
        assert(0 <= i && i < N);
        int k = i + offset;
        push_lazy_deep(k);
        dat[k] = ACT(f, dat[k]);
        pull_dat_deep(k);
    }
    // apply f for interval [l, r)
    void apply(int l, int r, const Action &f) {
        assert(0 <= l && l <= r && r <= N);
        if (l == r) return;
        l += offset, r += offset;
        for (int h = log; h >= 1; --h) {
            if (((l >> h) << h) != l) push_lazy(l >> h);
            if (((r >> h) << h) != r) push_lazy((r - 1) >> h);
        }
        int original_l = l, original_r = r;
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) apply_lazy(l++, f);
            if (r & 1) apply_lazy(--r, f);
        }
        l = original_l, r = original_r;
        for (int h = 1; h <= log; ++h) {
            if (((l >> h) << h) != l) pull_dat(l >> h);
            if (((r >> h) << h) != r) pull_dat((r - 1) >> h);
        }
    }
    
    // get prod of interval [l, r)
    Monoid prod(int l, int r) {
        assert(0 <= l && l <= r && r <= N);
        if (l == r) return IDENTITY_MONOID;
        l += offset, r += offset;
        for (int h = log; h >= 1; --h) {
            if (((l >> h) << h) != l) push_lazy(l >> h);
            if (((r >> h) << h) != r) push_lazy(r >> h);
        }
        Monoid val_left = IDENTITY_MONOID, val_right = IDENTITY_MONOID;
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) val_left = OP(val_left, dat[l++]);
            if (r & 1) val_right = OP(dat[--r], val_right);
        }
        return OP(val_left, val_right);
    }
    Monoid all_prod() {
        return dat[1];
    }
    
    // get max r that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int max_right(const function<bool(Monoid)> f, int l = 0) {
        if (l == N) return N;
        l += offset;
        push_lazy_deep(l);
        Monoid sum = IDENTITY_MONOID;
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(OP(sum, dat[l]))) {
                while (l < offset) {
                    push_lazy(l);
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
        push_lazy_deep(r - 1);
        Monoid sum = IDENTITY_MONOID;
        do {
            --r;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(OP(dat[r], sum))) {
                while (r < offset) {
                    push_lazy(r);
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
    
    // debug stream
    friend ostream& operator << (ostream &s, LazySegmentTree seg) {
        for (int i = 0; i < (int)seg.size(); ++i) {
            s << seg[i];
            if (i != (int)seg.size() - 1) s << " ";
        }
        return s;
    }
    
    // dump
    void dump() {
        for (int i = 0; i <= log; ++i) {
            for (int j = (1 << i); j < (1 << (i + 1)); ++j) {
                cout << "{" << dat[j] << "," << lazy[j] << "} ";
            }
            cout << endl;
        }
    }
};

// various segment trees
template<class Monoid> auto seg_op_min = [](Monoid p, Monoid q) { return min(p, q); };
template<class Monoid> auto seg_op_max = [](Monoid p, Monoid q) { return max(p, q); };
template<class Monoid> auto seg_op_sum = [](Monoid p, Monoid q) { return p + q; };
template<class Monoid> SegmentTree<Monoid> RangeMin(int N = 0) {
    return SegmentTree<Monoid>(N, seg_op_min<Monoid>, numeric_limits<Monoid>::max()/2);
}
template<class Monoid> SegmentTree<Monoid> RangeMin(const vector<Monoid> &v) {
    return SegmentTree<Monoid>(v, seg_op_min<Monoid>, numeric_limits<Monoid>::max()/2);
}
template<class Monoid> SegmentTree<Monoid> RangeMax(int N = 0) {
    return SegmentTree<Monoid>(N, seg_op_max<Monoid>, -numeric_limits<Monoid>::max()/2);
}
template<class Monoid> SegmentTree<Monoid> RangeMax(const vector<Monoid> &v) {
    return SegmentTree<Monoid>(v, seg_op_max<Monoid>, -numeric_limits<Monoid>::max()/2);
}
template<class Monoid> SegmentTree<Monoid> RangeSum(int N = 0) {
    return SegmentTree<Monoid>(N, seg_op_sum<Monoid>, Monoid(0));
}
template<class Monoid> SegmentTree<Monoid> RangeSum(const vector<Monoid> &v) {
    return SegmentTree<Monoid>(v, seg_op_sum<Monoid>, Monoid(0));
}

// various lazy segment trees
template<class Monoid> auto seg_op_sum_with_sum = [](pair<Monoid, long long> p, pair<Monoid, long long> q) { 
    return make_pair(p.first + q.first, p.second + q.second);
};
template<class Monoid, class Action> auto seg_act_change = [](Action f, Monoid x) {
    return (f != Action(-1) ? f : x);
};
template<class Monoid, class Action> auto seg_act_chmin = [](Action f, Monoid x) {
    return min(Monoid(f), x);
};
template<class Monoid, class Action> auto seg_act_chmax = [](Action f, Monoid x) {
    return max(Monoid(f), x);
};
template<class Monoid, class Action> auto seg_act_add = [](Action f, Monoid x) {
    return f + x;
};
template<class Monoid, class Action> auto seg_act_change_with_sum = [](Action f, pair<Monoid, long long> x) {
    return (f != Action(-1) ? make_pair(f * x.second, x.second) : x);
};
template<class Action> auto seg_comp_change = [](Action g, Action f) {
    return (g != Action(-1) ? g : f);
};
template<class Action> auto seg_comp_chmin = [](Action g, Action f) {
    return min(g, f);
};
template<class Action> auto seg_comp_chmax = [](Action g, Action f) {
    return max(g, f);
};
template<class Action> auto seg_comp_add = [](Action g, Action f) {
    return g + f;
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMin(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_min<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        numeric_limits<Monoid>::max()/2, Action(-1));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMin(const vector<Monoid> &v) {
    return LazySegmentTree<Monoid, Action>(
        v, seg_op_min<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        numeric_limits<Monoid>::max()/2, Action(-1));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMax(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_max<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        -numeric_limits<Monoid>::max()/2, Action(-1));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMax(const vector<Monoid> &v) {
    return LazySegmentTree<Monoid, Action>(
        v, seg_op_max<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        -numeric_limits<Monoid>::max()/2, Action(-1));
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid, long long>, Action> RangeChangeRangeSum(int N = 0) {
    return LazySegmentTree<pair<Monoid, long long>, Action>(
        N, seg_op_sum_with_sum<Monoid>, seg_act_change_with_sum<Monoid, Action>, seg_comp_change<Action>,
        make_pair(Monoid(0), 1), Action(-1));
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid, long long>, Action> RangeChangeRangeSum(const vector<Monoid> &v) {
    return LazySegmentTree<pair<Monoid, long long>, Action>(
        v, seg_op_sum_with_sum<Monoid>, seg_act_change_with_sum<Monoid, Action>, seg_comp_change<Action>,
        make_pair(Monoid(0), 1), Action(-1));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChminRangeMin(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_min<Monoid>, seg_act_chmin<Monoid, Action>, seg_comp_chmin<Action>,
        numeric_limits<Monoid>::max()/2, numeric_limits<Action>::max()/2);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChminRangeMin(const vector<Monoid> &v) {
    return LazySegmentTree<Monoid, Action>(
        v, seg_op_min<Monoid>, seg_act_chmin<Monoid, Action>, seg_comp_chmin<Action>,
        numeric_limits<Monoid>::max()/2, numeric_limits<Action>::max()/2);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChmaxRangeMax(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_max<Monoid>, seg_act_chmax<Monoid, Action>, seg_comp_chmax<Action>,
        -numeric_limits<Monoid>::max()/2, -numeric_limits<Action>::max()/2);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChmaxRangeMax(const vector<Monoid> &v) {
    return LazySegmentTree<Monoid, Action>(
        v, seg_op_max<Monoid>, seg_act_chmax<Monoid, Action>, seg_comp_chmax<Action>,
        -numeric_limits<Monoid>::max()/2, -numeric_limits<Action>::max()/2);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeAddRangeMin(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_min<Monoid>, seg_act_add<Monoid, Action>, seg_comp_add<Action>,
        numeric_limits<Monoid>::max()/2, Action(0));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeAddRangeMin(const vector<Monoid> &v) {
    return LazySegmentTree<Monoid, Action>(
        v, seg_op_min<Monoid>, seg_act_add<Monoid, Action>, seg_comp_add<Action>,
        numeric_limits<Monoid>::max()/2, Action(0));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeAddRangeMax(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_max<Monoid>, seg_act_add<Monoid, Action>, seg_comp_add<Action>,
        -numeric_limits<Monoid>::max()/2, Action(0));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeAddRangeMax(const vector<Monoid> &v) {
    return LazySegmentTree<Monoid, Action>(
        v, seg_op_max<Monoid>, seg_act_add<Monoid, Action>, seg_comp_add<Action>,
        -numeric_limits<Monoid>::max()/2, Action(0));
};


//------------------------------//
// Other Advanced Data Structures
//------------------------------//

// Bit Vector (for 64-bit non-negative integer)
struct BitVector {
    // block: bit vector
    // count: the number of 1 within each block
    unsigned int n, zeros;
    vector<unsigned long long> block;
    vector<unsigned int> count;
    
    // constructor
    BitVector() {}
    BitVector(const unsigned int num) {
        resize(num);
    }
    void resize(const unsigned int num) {
        n = num;
        block.assign(((num + 1) >> 6) + 1, 0);
        count.assign(block.size(), 0);
    }
    
    // set val(0 or 1) onto i-th bit, get i-th bit of val(0 or 1)
    void set(const unsigned int i, const unsigned long long val = 1LL) {
        assert((i >> 6) < block.size());
        block[i >> 6] |= (val << (i & 63));
    }
    unsigned int get(const unsigned int i) const {
        assert((i >> 6) < block.size());
        return (unsigned int)(block[i >> 6] >> (i & 63)) & 1;
    }
    void build() {
        for (unsigned int i = 1; i < block.size(); i++) {
            count[i] = count[i - 1] + __builtin_popcountll(block[i - 1]);
        }
        zeros = rank0(n);
    }
    
    // the number of 1 in [0, i)
    unsigned int rank1(const unsigned int i) const {
        assert((i >> 6) < count.size());
        assert((i >> 6) < block.size());
        return count[i >> 6] +
        __builtin_popcountll(block[i >> 6] & ((1ULL << (i & 63)) - 1ULL));
    }
    // the number of 1 in [i, j)
    unsigned int rank1(const unsigned int i, const unsigned int j) const {
        return rank1(j) - rank1(i);
    }
    // the number of 0 in [0, i)
    unsigned int rank0(const unsigned int i) const {
        return i - rank1(i);
    }
    // the number of 0 in [i, j)
    unsigned int rank0(const unsigned int i, const unsigned int j) const {
        return rank0(j) - rank0(i);
    }
    // the number of 0 in [0, n)
    unsigned int rank0() const {
        return zeros;
    }
};

// Wavelet Matrix (must vec[i] >= 0)
template<class T> struct WaveletMatrix {
    // inner data
    unsigned int n, height;
    vector<T> v;
    vector<BitVector> bv;
    vector<vector<long long>> sum;

    // constructor (sigma: the number of characters)
    WaveletMatrix() : n(0) {}
    WaveletMatrix(unsigned int n) : n(n), v(n) {}
    WaveletMatrix(const vector<T> &vec) : n(vec.size()), v(vec) {
        build();
    }
    void add(const T &val) {
        assert(v >= 0);
        v.push_back(v);
        n = v.size();
    }
    void set(unsigned int i, const T &val) {
        assert(i >= 0 && i < n && val >= 0);
        v[i] = val;
    }
    void build() {
        assert(n == (int)v.size());
        T mv = 1;
        for (int i = 0; i < n; ++i) mv = max(mv, v[i]);
        for (height = 1; mv != 0; mv >>= 1) ++height;
        vector<int> left(n), right(n), ord(n);
        iota(ord.begin(), ord.end(), 0);
        bv.assign(height, BitVector(n));
        sum.assign(height + 1, vector<long long>(n + 1, 0));
        for (int h = height - 1; h >= 0; --h) {
            int l = 0, r = 0;
            for (int i = 0; i < n; ++i) {
                if ((v[ord[i]] >> h) & 1) {
                    bv[h].set(i);
                    right[r++] = ord[i];
                } else {
                    left[l++] = ord[i];
                }
            }
            bv[h].build();
            ord.swap(left);
            for (int i = 0; i < r; ++i) ord[i + l] = right[i];
            for (int i = 0; i < n; ++i) sum[h][i + 1] = sum[h][i] + v[ord[i]];
        }
    }
    
    // access v[k]
    T access(int i) {
        T res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int i0 = bv[h].rank0(i);
            if (bv[h].get(i)) {
                i += bv[h].rank0() - i0;
                res |= T(1) << h;
            } else {
                i = i0;
            }
        }
        return res;
    }
    T operator [] (int i) {
        return access(i);
    }
    
    // count "i" s.t. v[i] = val, i \in [l, r)
    int rank(int l, int r, const T &val) {
        assert(0 <= l && l <= r && r <= n);
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if ((val >> h) & 1) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
            } else {
                l = l0;
                r = r0;
            }
        }
        return r - l;
    }
    
    // count "i" s.t. v[i] \in [lower, upper), i \in [l, r)
    int range_freq(int l, int r, const T &upper) {
        assert(0 <= l && l <= r && r <= n);
        int res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if ((upper >> h) & 1) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                res += r0 - l0;
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
    int range_freq(int l, int r, const T &lower, const T &upper) {
        return range_freq(l, r, upper) - range_freq(l, r, lower);
    }
    
    // the k-th (0-indexed) smallest value in [l, r)
    T k_th_smallest(int l, int r, int k) {
        assert(0 <= l && l <= r && r <= n);
        T res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if (r0 - l0 <= k) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                k -= r0 - l0;
                res |= T(1) << h;
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
    
    // the k-th (0-indexed) largest value in [l, r)
    T k_th_largest(int l, int r, int k) {
        assert(0 <= l && l <= r && r <= n);
       return k_th_smallest(l, r, r - l - k - 1);
    }
    
    // the sum of the top-k sum in [l, r)
    T top_k_sum(int l, int r, int k) {
        assert(0 <= l && l <= r && r <= n);
        if (l == r) return 0;
        T res = 0, val = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if (r0 - l0 <= k) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                k -= r0 - l0;
                val |= T(1) << h;
                res += sum[h][r0] - sum[h][l0];
            } else {
                l = l0;
                r = r0;
            }
        }
        res += val * k;
        return res;
    }
    
    // the max value (< val) in [l, r)
    T prev_value(int l, int r, T val) {
        assert(0 <= l && l <= r && r <= n);
        int num = range_freq(l, r, 0, val);
        if (num == 0) return T(-1);
        else return k_th_smallest(l, r, num - 1);
    }
    
    // the min value (>= val) in [l, r)
    T next_value(int l, int r, T val) {
        assert(0 <= l && l <= r && r <= n);
        int num = range_freq(l, r, 0, val);
        if (num == r - l) return T(-1);
        else return k_th_smallest(l, r, num);
    }
};

// 2D queries
template<class POS, class VAL> struct BITonWaveletMatrix {
    // inner data
    struct BIT {
        VAL UNITY_SUM = 0;
        int N;
        vector<VAL> dat;
        
        // [0, n)
        BIT() {}
        BIT(int n, VAL unity = 0) : UNITY_SUM(unity), N(n), dat(n, unity) { }
        void init(int n) {
            N = n;
            dat.assign(n, UNITY_SUM);
        }
        
        // a is 0-indexed
        void add(int a, VAL x) {
            for (int i = a; i < (int)dat.size(); i |= i + 1)
                dat[i] = dat[i] + x;
        }
        
        // [0, a), a is 0-indexed
        VAL sum(int a) {
            VAL res = UNITY_SUM;
            for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
                res = res + dat[i];
            return res;
        }
        
        // [a, b), a and b are 0-indexed
        VAL sum(int a, int b) {
            return sum(b) - sum(a);
        }
        
        // debug
        void print() {
            for (int i = 0; i < (int)dat.size(); ++i)
                cout << sum(i, i + 1) << ",";
            cout << endl;
        }
    };
    
    using Point = pair<POS, POS>;
    int n, height;
    POS mi = 0;
    vector<BitVector> bv;
    vector<Point> ps;
    vector<POS> ys;
    vector<BIT> bit;

    // constructor (sigma: the number of characters)
    // add_point() cannot be used after build()
    BITonWaveletMatrix() {}
    BITonWaveletMatrix(const vector<Point> &vec) {
        for (auto [x, y] : vec) add_point(x, y);
        build();
    }
    void add_point(POS x, POS y) {
        ps.emplace_back(x, y);
        ys.emplace_back(y);
        mi = min(mi, y);
    }
    int xid(POS x) const {
        return lower_bound(ps.begin(), ps.end(), Point(x, mi)) - ps.begin();
    }
    int yid(POS y) const {
        return lower_bound(ys.begin(), ys.end(), y) - ys.begin();
    }
    void build() {
        sort(ps.begin(), ps.end());
        ps.erase(unique(ps.begin(), ps.end()), ps.end());
        n = (int)ps.size();
        sort(ys.begin(), ys.end());
        ys.erase(unique(ys.begin(), ys.end()), ys.end());
        vector<int> v(n), left(n), right(n), ord(n);
        int mv = 1;
        for (int i = 0; i < n; ++i) {
            v[i] = yid(ps[i].second);
            mv = max(mv, v[i]);
        }
        for (height = 1; mv != 0; mv >>= 1) ++height;
        iota(ord.begin(), ord.end(), 0);
        bv.resize(height, BitVector(n));
        bit.assign(height + 1, BIT(n));
        for (int h = height - 1; h >= 0; --h) {
            int l = 0, r = 0;
            for (int i = 0; i < n; ++i) {
                if ((v[ord[i]] >> h) & 1) {
                    bv[h].set(i);
                    right[r++] = ord[i];
                } else {
                    left[l++] = ord[i];
                }
            }
            bv[h].build();
            ord.swap(left);
            for (int i = 0; i < r; ++i) ord[i + l] = right[i];
        }
    }
    
    // add
    void add(const POS x, const POS y, const VAL val) {
        int i = lower_bound(ps.begin(), ps.end(), Point(x, y)) - ps.begin();
        int j = yid(y);
        for (int h = height - 1; h >= 0; --h) {
            int i0 = bv[h].rank0(i);
            if ((j >> h) & 1) {
                i += bv[h].rank0() - i0;
            } else {
                i = i0;
            }
            bit[h].add(i, val);
        }
    }
    
    // sum
    VAL inner_sum(int l, int r, const POS upper) {
        assert(0 <= l && r <= n);
        VAL res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if ((upper >> h) & 1) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                res += bit[h].sum(l0, r0);
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
    VAL sum(const POS lx, const POS rx, const POS ly, const POS ry) {
        int l = xid(lx), r = xid(rx);
        return inner_sum(l, r, yid(ry)) - inner_sum(l, r, yid(ly));
    }
};


//------------------------------//
// String
//------------------------------//

// SA-IS (O(N))
template<class Str = string> struct SuffixArray {
    // data
    Str str;
    vector<int> sa;    // sa[i] : the starting index of the i-th smallest suffix (i = 0, 1, ..., n)
    vector<int> rank;  // rank[sa[i]] = i
    vector<int> lcp;   // lcp[i]: the lcp of sa[i] and sa[i+1] (i = 0, 1, ..., n-1)
    SparseTable<int> st;  // use for calcultating lcp(i, j)

    // getter
    int& operator [] (int i) { return sa[i]; }
    const int& operator [] (int i) const { return sa[i]; }
    vector<int> get_sa() { return sa; }
    vector<int> get_rank() { return rank; }
    vector<int> get_lcp() { return lcp; }

    // constructor
    SuffixArray() {}
    SuffixArray(const Str& str_, bool no_limit_elements = false) : str(str_) {
        build_sa(no_limit_elements);
    }
    void init(const Str& str_, bool no_limit_elements = false) {
        str = str_;
        build_sa(no_limit_elements);
    }
    void build_sa(bool no_limit_elements = false) {
        vector<int> s;
        int num_of_chars = 256;
        if (!no_limit_elements) {
            for (int i = 0; i < (int)str.size(); ++i) {
                s.push_back(str[i] + 1);
            }
        } else {
            unordered_map<int,int> dict;
            for (int i = 0; i < (int)str.size(); ++i) {
                if (!dict.count(str[i])) dict[str[i]] = dict.size();
            }
            for (int i = 0; i < (int)str.size(); ++i) {
                s.push_back(dict[str[i]] + 1);
            }
            num_of_chars = (int)dict.size();
        }
        s.push_back(0);
        sa = sa_is(s, num_of_chars);
        build_lcp(s);
        build_sparse_table();
    }

    // SA-IS
    // num_of_chars: # of characters
    vector<int> sa_is(vector<int> &s, int num_of_chars) {
        int N = (int)s.size();
        if (N == 0) return {};
        else if (N == 1) return {0};
        else if (N == 2) {
            if (s[0] < s[1]) return {0, 1};
            else return {1, 0};
        }

        vector<int> isa(N);
        vector<bool> ls(N, false);
        for (int i = N - 2; i >= 0; --i) {
            ls[i] = (s[i] == s[i + 1]) ? ls[i + 1] : (s[i] < s[i + 1]);
        }
        vector<int> sum_l(num_of_chars + 1, 0), sum_s(num_of_chars + 1, 0);
        for (int i = 0; i < N; ++i) {
            if (!ls[i]) ++sum_s[s[i]];
            else ++sum_l[s[i] + 1];
        }
        for (int i = 0; i <= num_of_chars; ++i) {
            sum_s[i] += sum_l[i];
            if (i < num_of_chars) sum_l[i + 1] += sum_s[i];
        }

        auto induce = [&](const vector<int> &lms) -> void {
            fill(isa.begin(), isa.end(), -1);
            vector<int> buf(num_of_chars + 1);
            copy(sum_s.begin(), sum_s.end(), buf.begin());
            for (auto d: lms) {
                if (d == N) continue;
                isa[buf[s[d]]++] = d;
            }
            copy(sum_l.begin(), sum_l.end(), buf.begin());
            isa[buf[s[N - 1]]++] = N - 1;
            for (int i = 0; i < N; ++i) {
                int v = isa[i];
                if (v >= 1 && !ls[v - 1]) {
                    isa[buf[s[v - 1]]++] = v - 1;
                }
            }
            copy(sum_l.begin(), sum_l.end(), buf.begin());
            for (int i = N - 1; i >= 0; --i) {
                int v = isa[i];
                if (v >= 1 && ls[v - 1]) {
                    isa[--buf[s[v - 1] + 1]] = v - 1;
                }
            }
        };
            
        vector<int> lms, lms_map(N + 1, -1);
        int M = 0;
        for (int i = 1; i < N; ++i) {
            if (!ls[i - 1] && ls[i]) {
                lms_map[i] = M++;
            }
        }
        lms.reserve(M);
        for (int i = 1; i < N; ++i) {
            if (!ls[i - 1] && ls[i]) {
                lms.push_back(i);
            }
        }
        induce(lms);

        if (M) {
            vector<int> lms2;
            lms2.reserve(isa.size());
            for (auto v: isa) {
                if (lms_map[v] != -1) lms2.push_back(v);
            }
            int rec_upper = 0;
            vector<int> rec_s(M);
            rec_s[lms_map[lms2[0]]] = 0;
            for (int i = 1; i < M; ++i) {
                int l = lms2[i - 1], r = lms2[i];
                int nl = (lms_map[l] + 1 < M) ? lms[lms_map[l] + 1] : N;
                int nr = (lms_map[r] + 1 < M) ? lms[lms_map[r] + 1] : N;
                bool same = true;
                if (nl - l != nr - r) same = false;
                else {
                    while (l < nl) {
                        if (s[l] != s[r]) break;
                        ++l, ++r;
                    }
                    if (l == N || s[l] != s[r]) same = false;
                }
                if (!same) ++rec_upper;
                rec_s[lms_map[lms2[i]]] = rec_upper;
            }
            auto rec_sa = sa_is(rec_s, rec_upper);

            vector<int> sorted_lms(M);
            for (int i = 0; i < M; ++i) {
                sorted_lms[i] = lms[rec_sa[i]];
            }
            induce(sorted_lms);
        }
        return isa;
    }

    // find min id that str.substr(sa[id]) >= T
    int lower_bound(const Str& T) {
        int left = -1, right = sa.size();
        while (right - left > 1) {
            int mid = (left + right) / 2;
            if (str.compare(sa[mid], string::npos, T) < 0)
                left = mid;
            else
                right = mid;
        }
        return right;
    }

    // find min id that str.substr(sa[id], T.size()) > T
    int upper_bound(const Str& T) {
        int left = -1, right = sa.size();
        while (right - left > 1) {
            int mid = (left + right) / 2;
            if (str.compare(sa[mid], T.size(), T) <= 0)
                left = mid;
            else
                right = mid;
        }
        return right;
    }
    
    // find min id that sa[id] >= str.substr(l, r-l)
    int lower_bound(int l, int r) {
        int left = -1, right = rank[l];
        while (right - left > 1) {
            int mid = (left + right) / 2;
            if (st.get(mid, rank[l]) < r - l) left = mid;
            else right = mid;
        }
        return right;
    }

    // search
    bool is_contain(const Str& T) {
        int lb = lower_bound(T);
        if (lb >= sa.size()) return false;
        return str.compare(sa[lb], T.size(), T) == 0;
    }

    // find lcp
    void build_lcp(const vector<int> &s) {
        int N = (int)s.size();
        rank.assign(N, 0), lcp.assign(N - 1, 0);
        for (int i = 0; i < N; ++i) rank[sa[i]] = i;
        int h = 0;
        for (int i = 0; i < N - 1; ++i) {
            int pi = sa[rank[i] - 1];
            if (h > 0) --h;
            for (; pi + h < N && i + h < N; ++h) {
                if (s[pi + h] != s[i + h]) break;
            }
            lcp[rank[i] - 1] = h;
        }
    }
    
    // build sparse table for calculating lcp
    void build_sparse_table() {
        st.init(lcp);
    }

    // calc lcp of str.sutstr(a) and str.substr(b)
    int get_lcp(int a, int b) {
        return st.get(min(rank[a], rank[b]), max(rank[a], rank[b]));
    }

    // debug
    void dump() {
        for (int i = 0; i < sa.size(); ++i) {
            cout << i << ": " << sa[i] << ", " << str.substr(sa[i]) << endl;
        }
    }
};

// Rolling Hash
template<class Str = string> struct RollingHash {
    static const int base1 = 1007, base2 = 2009;
    static const int mod1 = 1000000007, mod2 = 1000000009;
    vector<long long> hash1, hash2, power1, power2;

    // construct
    RollingHash(const Str &S) {
        int n = (int)S.size();
        hash1.assign(n+1, 0), hash2.assign(n+1, 0);
        power1.assign(n+1, 1), power2.assign(n+1, 1);
        for (int i = 0; i < n; ++i) {
            hash1[i+1] = (hash1[i] * base1 + S[i]) % mod1;
            hash2[i+1] = (hash2[i] * base2 + S[i]) % mod2;
            power1[i+1] = (power1[i] * base1) % mod1;
            power2[i+1] = (power2[i] * base2) % mod2;
        }
    }
    
    // get hash value of S[left:right]
    inline long long get(int l, int r) const {
        long long res1 = hash1[r] - hash1[l] * power1[r-l] % mod1;
        if (res1 < 0) res1 += mod1;
        long long res2 = hash2[r] - hash2[l] * power2[r-l] % mod2;
        if (res2 < 0) res2 += mod2;
        return res1 * mod2 + res2;
    }

    // get hash value of S
    inline long long get() const {
        return hash1.back() * mod2 + hash2.back();
    }

    // get lcp of S[a:] and S[b:]
    inline int getLCP(int a, int b) const {
        int len = min((int)hash1.size()-a, (int)hash1.size()-b);
        int low = 0, high = len;
        while (high - low > 1) {
            int mid = (low + high) >> 1;
            if (get(a, a+mid) != get(b, b+mid)) high = mid;
            else low = mid;
        }
        return low;
    }

    // get lcp of S[a:] and T[b:]
    inline int getLCP(const RollingHash &T, int a, int b) const {
        int len = min((int)hash1.size()-a, (int)hash1.size()-b);
        int low = 0, high = len;
        while (high - low > 1) {
            int mid = (low + high) >> 1;
            if (get(a, a+mid) != T.get(b, b+mid)) high = mid;
            else low = mid;
        }
        return low;
    }
};


//------------------------------//
// Geometry
//------------------------------//

// basic settings
long double EPS = 1e-10;  // to be set appropriately
constexpr long double PI = 3.141592653589793238462643383279502884L;
long double torad(long double deg) {return (long double)(deg) * PI / 180;}
long double todeg(long double ang) {return ang * 180 / PI;}

// Point or Vector
template<class DD> struct Point {
    // inner value
    DD x, y;
    
    // constructor
    constexpr Point() : x(0), y(0) {}
    constexpr Point(DD x, DD y) : x(x), y(y) {}
    
    // various functions
    constexpr Point conj() const {return Point(x, -y);}
    constexpr DD dot(const Point &r) const {return x * r.x + y * r.y;}
    constexpr DD cross(const Point &r) const {return x * r.y - y * r.x;}
    constexpr DD norm() const {return dot(*this);}
    constexpr long double abs() const {return sqrt(norm());}
    constexpr long double arg() const {
        if (this->eq(Point(0, 0))) return 0L;
        else if (x < -EPS && this->eq(Point(x, 0))) return PI;
        else return atan2((long double)(y), (long double)(x));
    }
    constexpr bool eq(const Point &r) const {return (*this - r).abs() <= EPS;}
    constexpr int sign() const {
        if (x >= -EPS && x <= EPS && y >= -EPS && y <= EPS) return 0;
        else if (y < -EPS || (y >= -EPS && y <= EPS && x > EPS)) return -1;
        else return 1;
    }
    constexpr Point rot90() const {return Point(-y, x);}
    constexpr Point rot(long double ang) const {
        return Point(cos(ang) * x - sin(ang) * y, sin(ang) * x + cos(ang) * y);
    }
    
    // arithmetic operators
    constexpr Point operator - () const {return Point(-x, -y);}
    constexpr Point operator + (const Point &r) const {return Point(*this) += r;}
    constexpr Point operator - (const Point &r) const {return Point(*this) -= r;}
    constexpr Point operator * (const Point &r) const {return Point(*this) *= r;}
    constexpr Point operator / (const Point &r) const {return Point(*this) /= r;}
    constexpr Point operator * (DD r) const {return Point(*this) *= r;}
    constexpr Point operator / (DD r) const {return Point(*this) /= r;}
    constexpr Point& operator += (const Point &r) {
        x += r.x, y += r.y;
        return *this;
    }
    constexpr Point& operator -= (const Point &r) {
        x -= r.x, y -= r.y;
        return *this;
    }
    constexpr Point& operator *= (const Point &r) {
        DD tx = x, ty = y;
        x = tx * r.x - ty * r.y;
        y = tx * r.y + ty * r.x;
        return *this;
    }
    constexpr Point& operator /= (const Point &r) {
        return *this *= r.conj() / r.norm();
    }
    constexpr Point& operator *= (DD r) {
        x *= r, y *= r;
        return *this;
    }
    constexpr Point& operator /= (DD r) {
        x /= r, y /= r;
        return *this;
    }

    // friend functions
    friend ostream& operator << (ostream &s, const Point &p) {
        return s << '(' << p.x << ", " << p.y << ')';
    }
    friend constexpr Point conj(const Point &p) {return p.conj();}
    friend constexpr DD dot(const Point &p, const Point &q) {return p.dot(q);}
    friend constexpr DD cross(const Point &p, const Point &q) {return p.cross(q);}
    friend constexpr DD norm(const Point &p) {return p.norm();}
    friend constexpr long double abs(const Point &p) {return p.abs();}
    friend constexpr long double arg(const Point &p) {return p.arg();}
    friend constexpr bool eq(const Point &p, const Point &q) {return p.eq(q);}
    friend constexpr int sign(const Point &p) {return p.sign();}
    friend constexpr Point rot90(const Point &p) {return p.rot90();}
    friend constexpr Point rot(const Point &p, long long ang) {return p.rot(ang);}
};

// necessary for some functions
template<class DD> constexpr bool operator < (const Point<DD> &p, const Point<DD> &q) {
    return (abs(p.x - q.x) > EPS ? p.x < q.x : p.y < q.y);
}

// Line
template<class DD> struct Line : vector<Point<DD>> {
    Line(Point<DD> a = Point<DD>(0, 0), Point<DD> b = Point<DD>(0, 0)) {
        this->push_back(a);
        this->push_back(b);
    }
    friend ostream& operator << (ostream &s, const Line<DD> &l) {
        return s << '{' << l[0] << ", " << l[1] << '}';
    }
};

// Circle
template<class DD> struct Circle : Point<DD> {
    DD r;
    Circle(Point<DD> p = Point<DD>(0, 0), DD r = 0) : Point<DD>(p), r(r) {}
    friend ostream& operator << (ostream &s, const Circle<DD> &c) {
        return s << '(' << c.x << ", " << c.y << ", " << c.r << ')';
    }
};

// arg sort
// by defining comparison
template<class DD> bool arg_cmp(const Point<DD> &p, const Point<DD> &q) {
    if (sign(p) != sign(q)) return sign(p) < sign(q);
    return (abs(cross(p, q)) > EPS ? cross(p, q) > EPS : norm(p) < norm(q));
}
template<class DD> void arg_sort(vector<Point<DD>> &v) {
    sort(v.begin(), v.end(), [&](const Point<DD> &p, const Point<DD> &q){return arg_cmp(p, q);});
}


//------------------------------//
// Tree
//------------------------------//

// Run Tree (including Euler Tour)
template<class T = long long> struct RunTree {
    // id[v][w] := the index of node w in G[v]
    vector<unordered_map<int, int>> id;

    // siz[v][i] := the size of subtree of G[v][i] with parent v
    vector<vector<long long>> siz;
    
    // for finding lca
    int root;
    vector<vector<int>> parent;
    vector<int> depth;

    // Euler tour
    vector<int> tour; // the node-number of i-th element of Euler-tour
    vector<int> v_s_id, v_t_id; // the index of Euler-tour of node v
    vector<int> e_id; // the index of edge e (v*2 + (0: root to leaf, 1: leaf to root))

    // constructor
    RunTree() {}
    RunTree(const Graph<T> &G, int root = 0) : root(root) {
        init(G, root);
    }
    
    // init
    void init(const Graph<T> &G, int root = 0) {
        int N = (int)G.size();
        id.assign(N, unordered_map<int,int>()), siz.assign(N, vector<long long>());
        for (int v = 0; v < N; v++) siz[v].assign((int)G[v].size(), 0);
        int h = 1, ord = 0;
        while ((1<<h) < N) h++;
        parent.assign(h, vector<int>(N, -1)), depth.resize(N);
        tour.resize(N*2-1), v_s_id.resize(N), v_t_id.resize(N), e_id.resize(N*2);
        rec(G, root, -1, 0, ord);
        for (int i = 0; i+1 < (int)parent.size(); ++i) {
            for (int v = 0; v < N; v++)
                if (parent[i][v] != -1)
                    parent[i+1][v] = parent[i][parent[i][v]];
        }
    }

    // get_size(u, v) := the size of subtree v with parent u
    long long get_size(int u, int v) {
        return siz[u][id[u][v]];
    }

    // get first / last id of node v in Euler tour
    int vs(int v) { return v_s_id[v]; }
    int vt(int v) { return v_t_id[v]; }
    int get_v(int id) { return tour[id]; }

    // get edge-id of (pv, v) in Euler tour
    int e(int v, bool leaf_to_root = false) {
        assert(v != root);
        if (!leaf_to_root) return e_id[v * 2];
        else return e_id[v * 2 + 1];
    }
    int e(int u, int v) {
        if (depth[u] < depth[v]) return e(v);
        else return e(u, false);
    }
    pair<int, int> get_e(int id) { 
        return make_pair(tour[id], tour[id + 1]);
    }

    // get_parent(v, p) := the parent of v directed for p
    int get_parent(int v) { return parent[0][v]; }
    int kth_ancestor(int v, int k) {
        if (k > depth[v]) return root;
        int goal_depth = depth[v] - k;
        for (int i = (int)parent.size()-1; i >= 0; i--)
            if (parent[i][v] != -1 && depth[parent[i][v]] >= goal_depth) 
                v = parent[i][v];
        return v;
    }
    int get_parent(int v, int p) {
        if (v == p) return -1;
        int lca = get_lca(v, p);
        if (lca != v) return parent[0][v];
        else return kth_ancestor(p, depth[p] - depth[v] - 1);
    }

    // lca(u, v)
    int get_lca(int u, int v) {
        if (depth[u] > depth[v]) swap(u, v);
        for (int i = 0; i < (int)parent.size(); i++) {
            if ((depth[v] - depth[u]) & (1<<i))
                v = parent[i][v];
        }
        if (u == v) return u;
        for (int i = (int)parent.size()-1; i >= 0; i--) {
            if (parent[i][u] != parent[i][v]) {
                u = parent[i][u];
                v = parent[i][v];
            }
        }
        return parent[0][u];
    }

    // dist(u, v)
    long long get_dist(int u, int v) {
        int lca = get_lca(u, v);
        return depth[u] + depth[v] - depth[lca]*2;
    }

    // is node v in s-t path?
    bool is_on_path(int s, int t, int v) {
        return get_dist(s, v) + get_dist(v, t) == get_dist(s, t);
    };
    
    // rec
    int rec(const Graph<T> &G, int v, int p, int d, int &ord) {
        int p_index = -1;
        int sum = 1;
        parent[0][v] = p, depth[v] = d;
        tour[ord] = v, v_s_id[v] = v_t_id[v] = ord;
        ord++;
        for (int i = 0; i < (int)G[v].size(); i++) {
            int ch = G[v][i].to;
            id[v][ch] = i;
            if (ch == p) {
                p_index = i;
                continue;
            }
            e_id[ch * 2] = ord - 1;
            int s = rec(G, ch, v, d+1, ord);
            siz[v][i] = s;
            sum += s;
            tour[ord] = v;
            v_t_id[v] = ord;
            e_id[ch * 2 + 1] = ord - 1;
            ord++;
        }
        if (p_index != -1) siz[v][p_index] = (int)G.size() - sum;
        return sum;
    }
};

// find diameter of tree
template<class Weight = long long> pair<Weight, vector<Edge<Weight>>> calc_diameter
(const Graph<Weight> &G) {
    Weight length = 0;
    vector<Edge<Weight>> path;
    vector<Edge<Weight>> prev(G.size(), Edge<Weight>(-1, -1));

    auto dfs = [&](auto &&dfs, int v, int p, bool record = true) -> pair<int, Weight> {
        pair<int, Weight> res{v, 0};
        for (const auto &e : G[v]) {
            if (e.to == p) continue;
            auto tmp = dfs(dfs, e.to, v, record);
            tmp.second += e.val;
            if (tmp.second > res.second) {
                res = tmp;
                if (record) prev[e.to] = e;
            }
        }
        return res;
    };

    auto [leaf, distance] = dfs(dfs, 0, -1, false);
    prev.assign((int)G.size(), Edge<Weight>(-1, -1));
    auto [most_distant_v, distance2] = dfs(dfs, leaf, -1, true);
    int cur = most_distant_v;
    while (cur != -1) {
        const auto &e = prev[cur];
        if (e.from == -1) break;
        length += e.val, path.emplace_back(e);
        cur = e.from;
    }
    reverse(path.begin(), path.end());
    return {length, path};
}


//------------------------------//
// Solver
//------------------------------//

void output_bigcase(string filepath = "./big.in") {
	freopen( filepath.c_str(), "w", stdout );

    // テストケース生成
}

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    
    
}








