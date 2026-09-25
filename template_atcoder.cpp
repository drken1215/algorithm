// code template is in https://github.com/drken1215/algorithm/blob/master/template_atcoder.cpp
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

// Yes/No
void yes(bool a) { cout << (a ? "yes" : "no") << endl; }
void YES(bool a) { cout << (a ? "YES" : "NO") << endl; }
void Yes(bool a) { cout << (a ? "Yes" : "No") << endl; }

// neighbors
const vector<int> DX = {1, 0, -1, 0, 1, -1, 1, -1};
const vector<int> DY = {0, 1, 0, -1, 1, -1, -1, 1};

// Edge Class
template<class T = long long> struct Edge {
    int from, to;
    T val;
    Edge() : from(-1), to(-1) { }
    Edge(int f, int t, T v = 1) : from(f), to(t), val(v) {}
    Edge(const Edge&) = default;
    Edge& operator = (const Edge&) = default;
    friend ostream& operator << (ostream& s, const Edge& e) {
        return s << e.from << "->" << e.to << "(" << e.val << ")";
    }
};

// graph class
template<class T = long long> struct Graph {
    int V;
    bool record_reversed_edges = false, record_edge_index = false;
    vector<vector<Edge<T>>> list;
    vector<vector<Edge<T>>> reversed_list;
    vector<unordered_map<int, int>> id;  // id[v][w] := the index of node w in G[v]

    // constructors
    Graph(int n = 0, bool rre = false, bool rei = false) {
        init(n, rre, rei);
    }
    void init(int n = 0, bool rre = false, bool rei = false) {
        V = n, record_reversed_edges = rre, record_edge_index = rei;
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
        if (record_reversed_edges) reversed_list[to].push_back(Edge(to, from, val));
        if (from != to) {
            if (record_edge_index) id[to][from] = (int)list[to].size(); 
            list[to].push_back(Edge(to, from, val));
            if (record_reversed_edges) reversed_list[from].push_back(Edge(from, to, val));
        }
    }

    // input (only tree-case)
    friend istream& operator >> (istream &is, Graph &G) {
        for (int i = 0; i < G.V - 1; i++) {
            int u, v;
            is >> u >> v, u--, v--;
            G.add_bidirected_edge(u, v);
        }
        return is;
    }

    // output
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
template<typename T, typename... Ts>
vector<T> merge_and_unique(const vector<T> &a, const Ts &...z) {
    vector<T> res = a;
    (res.insert(res.end(), z.begin(), z.end()), ...);
    sort(res.begin(), res.end());
    res.erase(unique(res.begin(), res.end()), res.end());
    return res;
}

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
// mod algorithms
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
const int CHANGE_UNITY = 123456789;
template<class Monoid> auto seg_op_sum_with_sum = [](pair<Monoid, long long> p, pair<Monoid, long long> q) { 
    return make_pair(p.first + q.first, p.second + q.second);
};
template<class Monoid, class Action> auto seg_act_change = [](Action f, Monoid x) {
    return (f != Action(CHANGE_UNITY) ? f : x);
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
    return (f != Action(CHANGE_UNITY) ? make_pair(f * x.second, x.second) : x);
};
template<class Action> auto seg_comp_change = [](Action g, Action f) {
    return (g != Action(CHANGE_UNITY) ? g : f);
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
        numeric_limits<Monoid>::max()/2, Action(CHANGE_UNITY));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMin(const vector<Monoid> &v) {
    return LazySegmentTree<Monoid, Action>(
        v, seg_op_min<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        numeric_limits<Monoid>::max()/2, Action(CHANGE_UNITY));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMax(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_max<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        -numeric_limits<Monoid>::max()/2, Action(CHANGE_UNITY));
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMax(const vector<Monoid> &v) {
    return LazySegmentTree<Monoid, Action>(
        v, seg_op_max<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        -numeric_limits<Monoid>::max()/2, Action(CHANGE_UNITY));
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid, long long>, Action> RangeChangeRangeSum(int N = 0) {
    return LazySegmentTree<pair<Monoid, long long>, Action>(
        N, seg_op_sum_with_sum<Monoid>, seg_act_change_with_sum<Monoid, Action>, seg_comp_change<Action>,
        make_pair(Monoid(0), 0), Action(CHANGE_UNITY));
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid, long long>, Action> RangeChangeRangeSum(const vector<Monoid> &v) {
    vector<pair<Monoid, long long>> ini(v.size());
    for (int i = 0; i < (int)v.size(); i++) ini[i] = make_pair(v[i], 1);
    return LazySegmentTree<pair<Monoid, long long>, Action>(
        ini, seg_op_sum_with_sum<Monoid>, seg_act_change_with_sum<Monoid, Action>, seg_comp_change<Action>,
        make_pair(Monoid(0), 0), Action(CHANGE_UNITY));
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








