// code template is in https://github.com/drken1215/algorithm/blob/master/template_contest.cpp
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;

template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }

using pint = pair<int, int>;
using pll = pair<long long, long long>;
using tint = array<int, 3>;
using tll = array<long long, 3>;
using fint = array<int, 4>;
using fll = array<long long, 4>;
using qint = array<int, 5>;
using qll = array<long long, 5>;
using vint = vector<int>;
using vll = vector<long long>;
using ll = long long;
using u32 = unsigned int;
using u64 = unsigned long long;
using i128 = __int128_t;
using u128 = __uint128_t;
template <class T>
using min_priority_queue = priority_queue<T, vector<T>, greater<T>>;

#define REP(i, a) for (long long i = 0; i < (long long)(a); i++)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); i++)
#define RREP(i, a) for (long long i = (a)-1; i >= (long long)(0); --i)
#define RREP2(i, a, b) for (long long i = (b)-1; i >= (long long)(a); --i)
#define EB emplace_back
#define PB push_back
#define MP make_pair
#define MT make_tuple
#define FI first
#define SE second
#define ALL(x) x.begin(), x.end()
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl

// debug stream
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, array<T, 3> P)
{ return s << '<' << P[0] << ", " << P[1] << ", " << P[2] << '>'; }
template<class T> ostream& operator << (ostream &s, array<T, 4> P)
{ return s << '<' << P[0] << ", " << P[1] << ", " << P[2] << ", " << P[3] << '>'; }
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, deque<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, vector<vector<T> > P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, set<T> P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, multiset<T> P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, unordered_set<T> P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T1, class T2> ostream& operator << (ostream &s, map<T1,T2> P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
template<class T1, class T2> ostream& operator << (ostream &s, unordered_map<T1,T2> P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }


/*///////////////////////////////////////////////////////*/
// Fast IO
/*///////////////////////////////////////////////////////*/

struct FastRead {
    static constexpr int BUF_SIZE = 1 << 17;

private:
    FILE *stream_;
    array<char, BUF_SIZE> buf_;
    char *begin_, *end_, *ptr_;

    // reader
    void skip_space() {
        while (*ptr_ <= ' ') ++ptr_;
    }
    template<int N = 0> void read() {
        if (const auto n = end_ - ptr_; n <= N) {
            ignore = fread(copy_n(ptr_, n, begin_), 1, BUF_SIZE - n, stream_);
            ptr_ = begin_;
        }
    }
    
    // parser
    template<typename T> void parse(T &x) {
        common_type_t<T, uint64_t> x2 = 0;
        while (true) {
            uint64_t v;
            memcpy(&v, ptr_, 8);
            if ((v -= 0x3030303030303030) & 0x8080808080808080) break;
            v = (v * 10 + (v >> 8)) & 0xff00ff00ff00ff;
            v = (v * 100 + (v >> 16)) & 0xffff0000ffff;
            v = (v * 10000 + (v >> 32)) & 0xffffffff;
            x2 = 100000000 * x2 + v;
            ptr_ += 8;
        }
        while (true) {
            uint32_t v;
            memcpy(&v, ptr_, 4);
            if ((v -= 0x30303030) & 0x80808080) break;
            v = (v * 10 + (v >> 8)) & 0xff00ff;
            v = (v * 100 + (v >> 16)) & 0xffff;
            x2 = 10000 * x2 + v;
            ptr_ += 4;
            break;
        }
        while (true) {
            uint16_t v;
            memcpy(&v, ptr_, 2);
            if ((v -= 0x3030) & 0x8080) break;
            v = (v * 10 + (v >> 8)) & 0xff;
            x2 = 100 * x2 + v;
            ptr_ += 2;
            break;
        }
        if (' ' < *ptr_) {
            x2 *= 10;
            x2 += *ptr_++ - '0';
        }
        ++ptr_;
        x = static_cast<T>(x2);
    }
    
public:
    // constructor
    FastRead() : FastRead(stdin) {}
    explicit FastRead(const filesystem::path& p) : FastRead(fopen(p.c_str(), "r")) {}
    explicit FastRead(FILE *stream)
    : stream_(stream), begin_(buf_.data()), end_(begin_ + BUF_SIZE), ptr_(end_) { 
        read(); 
    }
    ~FastRead() { 
        if (stream_ != stdin) fclose(stream_); 
    }
    FastRead(const FastRead&) = delete;
    FastRead &operator = (const FastRead&) = delete;
    
    // operators
    template<unsigned_integral T> void operator () (T &x) {
        skip_space();
        read<64>();
        parse(x);
    }
    template<signed_integral T> void operator () (T &x) {
        skip_space();
        read<64>();
        make_unsigned_t<T> u;
        if (*ptr_ == '-') {
            ++ptr_;
            parse(u);
            u = -u;
        } else {
            parse(u);
        }
        x = u;
    }
    void operator () (char &x) {
        skip_space();
        read<64>();
        x = *ptr_;
        ++ptr_;
    }
    void operator () (string &x) {
        x = "";
        skip_space();
        read<64>();
        while (*ptr_ > ' ' && *ptr_ != '\0') {
            x.push_back(*ptr_);
            ++ptr_;
        }
        ++ptr_;
    }
    template<class... Ts> requires(sizeof...(Ts) != 1) void operator () (Ts&... xs) {
        ((*this)(xs), ...);
    }
    template<class T> FastRead& operator >> (T &x) { (*this)(x); return *this; }
};

class FastWrite {
    static constexpr int BUF_SIZE = 1 << 17;

private:
    FILE *stream_;
    array<char, BUF_SIZE> buf_;
    char *begin_, *end_, *ptr_;
    
    // preparation
    template <class T> static constexpr int DIGITS = numeric_limits<T>::digits10 + 1;
    template <class T> static constexpr auto POW10 = [] {
        array<T, DIGITS<T>> ret;
        ret[0] = 1;
        for (int i = 1; i < DIGITS<T>; ++i) {
            ret[i] = 10 * ret[i - 1];
        }
        return ret;
    } ();
    static constexpr auto LUT = [] {
        array<char, 40000> res;
        char* p = res.data();
        char a = '0', b = '0', c = '0', d = '0';
        do {
            *p++ = a, *p++ = b, *p++ = c, *p++ = d;
        } while (d++ < '9'
                 || (d = '0', c++ < '9'
                     || (c = '0', b++ < '9'
                         || (b = '0', a++ < '9'))));
        return res;
    } ();
    
    // flush
    template<int N = BUF_SIZE> void flush() {
        if (end_ - ptr_ <= N) {
            fwrite(begin_, 1, ptr_ - begin_, stream_);
            ptr_ = begin_;
        }
    }
    
    // writer
    template<int N = 4> void le4(uint64_t x) {
        if constexpr (1 < N) {
            if (x < POW10<uint64_t>[N - 1]) {
                le4<N - 1>(x);
                return;
            }
        }
        ptr_ = copy_n(&LUT[x * 4 + (4 - N)], N, ptr_);
    }
    template<int N> void w4(uint64_t x) {
        if constexpr (0 < N) {
            ptr_ = copy_n(&LUT[x / POW10<uint64_t>[N - 4] * 4], 4, ptr_);
            w4<N - 4>(x % POW10<uint64_t>[N - 4]);
        }
    }
    template<int N> void write(uint64_t x) {
        if constexpr (N < DIGITS<uint64_t>) {
            if (POW10<uint64_t>[N] <= x) {
                write<N + 4>(x);
                return;
            }
        }
        le4(x / POW10<uint64_t>[N - 4]);
        w4<N - 4>(x % POW10<uint64_t>[N - 4]);
    }
    template<typename T> void write(T x) {
        write<4>(x);
    }
    void write(__uint128_t x) {
        if (x < POW10<__uint128_t>[16]) {
            write(static_cast<uint64_t>(x));
        } else if (x < POW10<__uint128_t>[32]) {
            write(static_cast<uint64_t>(x / POW10<__uint128_t>[16]));
            w4<16>(static_cast<uint64_t>(x % POW10<__uint128_t>[16]));
        } else {
            write(static_cast<uint64_t>(x / POW10<__uint128_t>[32]));
            x %= POW10<__uint128_t>[32];
            w4<16>(static_cast<uint64_t>(x / POW10<__uint128_t>[16]));
            w4<16>(static_cast<uint64_t>(x % POW10<__uint128_t>[16]));
        }
    }
    
public:
    // constructor
    FastWrite() : FastWrite(stdout) {}
    explicit FastWrite(const filesystem::path& p) : FastWrite(fopen(p.c_str(), "w")) {}
    explicit FastWrite(FILE* stream)
    : stream_(stream), begin_(buf_.data()), end_(begin_ + BUF_SIZE), ptr_(begin_) {}
    ~FastWrite() {
        flush();
        if (stream_ != stdout) { fclose(stream_); }
    }
    FastWrite(const FastWrite&) = delete;
    FastWrite& operator = (const FastWrite&) = delete;
    
    // operators
    template<unsigned_integral T> void operator () (T x) {
        flush<DIGITS<T>>();
        write(x);
    }
    template<signed_integral T> void operator () (T x) {
        flush<1 + DIGITS<T>>();
        using U = make_unsigned_t<T>;
        const U u = x;
        if (x < 0) {
            *ptr_++ = '-';
            write(static_cast<U>(-u));
        } else {
            write(u);
        }
    }
    void operator () (char c) {
        flush<1>();
        *ptr_++ = c;
    }
    void operator () (string_view s) {
        while (!s.empty()) {
            flush<0>();
            const auto n = min(ssize(s), end_ - ptr_);
            if (n == BUF_SIZE) {
                fwrite(s.data(), 1, BUF_SIZE, stream_);
            } else {
                ptr_ = copy_n(s.data(), n, ptr_);
            }
            s.remove_prefix(n);
        }
        flush<0>();
    }
    template <char End = '\n', char Sep = ' ', class T, class... Ts>
    void ln(T&& x, Ts&&... xs) {
        (*this)(std::forward<T>(x));
        if constexpr (sizeof...(Ts) == 0) {
            *ptr_++ = End;
        } else {
            *ptr_++ = Sep;
            ln<End, Sep>(std::forward<Ts>(xs)...);
        }
    }
    template<class T> FastWrite& operator << (T x) { (*this)(x); return *this; }
};


/*///////////////////////////////////////////////////////*/
// Utility
/*///////////////////////////////////////////////////////*/

// 4-neighbor
const vector<int> dx = {1, 0, -1, 0};
const vector<int> dy = {0, 1, 0, -1};

// 8-neighbor
const vector<int> dx8 = {1, 0, -1, 0, 1, -1, 1, -1};
const vector<int> dy8 = {0, 1, 0, -1, 1, 1, -1, -1};

// bit
int popcnt(int x) { return __builtin_popcount(x); }
int popcnt(u32 x) { return __builtin_popcount(x); }
int popcnt(ll x) { return __builtin_popcountll(x); }
int popcnt(u64 x) { return __builtin_popcountll(x); }

// xor128による乱数生成、周期は2^128-1
unsigned int randInt() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty; ty = tz; tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}
int randInt(int minv, int maxv) {
    return randInt() % (maxv - minv + 1) + minv;
}
ll randll(ll minv, ll maxv) {
    ll a = randInt(), b = randInt();
    return (a * (1LL<<29) + b) % (maxv - minv + 1) + minv;
}
template<class T> void shuffle(vector<T>& vec) {
    int n = vec.size();
    for (int i = n - 1; i > 0; --i) {
        int k = randInt() % (i + 1);
        swap(vec[i], vec[k]);
    }
}

// floor, ceil
template<class T> T floor(T a, T b) {
    if (a % b == 0 || a >= 0) return a / b;
    else return -((-a) / b) - 1;
}
template<class T> T ceil(T x, T y) {
    return floor(x + y - 1, y);
}

// kth root
// N < 2^64, K <= 64
uint64_t kth_root(uint64_t N, uint64_t K) {
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

// int 128
constexpr i128 to_integer(const string &s) {
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


/*/////////////////////////////*/
// eratostenes
/*/////////////////////////////*/

// isprime[n] := is n prime?
// mebius[n] := mebius value of n
// min_factor[n] := the min prime-factor of n
// euler[n] := euler function value of n
struct Eratos {
    vector<int> primes;
    vector<bool> isprime;
    vector<int> mebius, min_factor, euler;

    // constructor, getter
    Eratos(int MAX) : primes(),
                      isprime(MAX+1, true),
                      mebius(MAX+1, 1),
                      min_factor(MAX+1, -1),
                      euler(MAX+1) {
        isprime[0] = isprime[1] = false;
        min_factor[0] = 0, min_factor[1] = 1;
        for (int i = 1; i <= MAX; i++) euler[i] = i;
        for (int i = 2; i <= MAX; ++i) {
            if (!isprime[i]) continue;
            primes.push_back(i);
            mebius[i] = -1;
            min_factor[i] = i;
            euler[i] = i - 1;
            for (int j = i*2; j <= MAX; j += i) {
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
        vector<pair<int,int> > res;
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


/*/////////////////////////////*/
// modint, FPS
/*/////////////////////////////*/

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

// dynamic modint
struct DynamicModint {
    using mint = DynamicModint;
    
    // static menber
    static int MOD;
    
    // inner value
    long long val;
    
    // constructor
    DynamicModint() : val(0) { }
    DynamicModint(long long v) : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    long long get() const { return val; }
    static int get_mod() { return MOD; }
    static void set_mod(int mod) { MOD = mod; }
    
    // arithmetic operators
    mint operator + () const { return mint(*this); }
    mint operator - () const { return mint(0) - mint(*this); }
    mint operator + (const mint &r) const { return mint(*this) += r; }
    mint operator - (const mint &r) const { return mint(*this) -= r; }
    mint operator * (const mint &r) const { return mint(*this) *= r; }
    mint operator / (const mint &r) const { return mint(*this) /= r; }
    mint& operator += (const mint &r) {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    mint& operator -= (const mint &r) {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    mint& operator *= (const mint &r) {
        val = val * r.val % MOD;
        return *this;
    }
    mint& operator /= (const mint &r) {
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
    mint pow(long long n) const {
        mint res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    mint inv() const {
        mint res(1), div(*this);
        return res / div;
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
        is >> x.val;
        x.val %= x.get_mod();
        if (x.val < 0) x.val += x.get_mod();
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

// FPS
namespace NTT {
    long long modpow(long long a, long long n, int mod) {
        long long res = 1;
        while (n > 0) {
            if (n & 1) res = res * a % mod;
            a = a * a % mod;
            n >>= 1;
        }
        return res;
    }

    long long modinv(long long a, int mod) {
        long long b = mod, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        u %= mod;
        if (u < 0) u += mod;
        return u;
    }

    int calc_primitive_root(int mod) {
        if (mod == 2) return 1;
        if (mod == 167772161) return 3;
        if (mod == 469762049) return 3;
        if (mod == 754974721) return 11;
        if (mod == 998244353) return 3;
        int divs[20] = {};
        divs[0] = 2;
        int cnt = 1;
        long long x = (mod - 1) / 2;
        while (x % 2 == 0) x /= 2;
        for (long long i = 3; i * i <= x; i += 2) {
            if (x % i == 0) {
                divs[cnt++] = i;
                while (x % i == 0) x /= i;
            }
        }
        if (x > 1) divs[cnt++] = x;
        for (int g = 2;; g++) {
            bool ok = true;
            for (int i = 0; i < cnt; i++) {
                if (modpow(g, (mod - 1) / divs[i], mod) == 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) return g;
        }
    }

    int get_fft_size(int N, int M) {
        int size_a = 1, size_b = 1;
        while (size_a < N) size_a <<= 1;
        while (size_b < M) size_b <<= 1;
        return max(size_a, size_b) << 1;
    }

    // number-theoretic transform
    template<class mint> void trans(vector<mint> &v, bool inv = false) {
        if (v.empty()) return;
        int N = (int)v.size();
        int MOD = v[0].get_mod();
        int PR = calc_primitive_root(MOD);
        static bool first = true;
        static vector<long long> vbw(30), vibw(30);
        if (first) {
            first = false;
            for (int k = 0; k < 30; ++k) {
                vbw[k] = modpow(PR, (MOD - 1) >> (k + 1), MOD);
                vibw[k] = modinv(vbw[k], MOD);
            }
        }
        for (int i = 0, j = 1; j < N - 1; j++) {
            for (int k = N >> 1; k > (i ^= k); k >>= 1);
            if (i > j) swap(v[i], v[j]);
        }
        for (int k = 0, t = 2; t <= N; ++k, t <<= 1) {
            long long bw = vbw[k];
            if (inv) bw = vibw[k];
            for (int i = 0; i < N; i += t) {
                mint w = 1;
                for (int j = 0; j < t/2; ++j) {
                    int j1 = i + j, j2 = i + j + t/2;
                    mint c1 = v[j1], c2 = v[j2] * w;
                    v[j1] = c1 + c2;
                    v[j2] = c1 - c2;
                    w *= bw;
                }
            }
        }
        if (inv) {
            long long invN = modinv(N, MOD);
            for (int i = 0; i < N; ++i) v[i] = v[i] * invN;
        }
    }

    // for garner
    static constexpr int MOD0 = 754974721;
    static constexpr int MOD1 = 167772161;
    static constexpr int MOD2 = 469762049;
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    static const mint1 imod0 = 95869806; // modinv(MOD0, MOD1);
    static const mint2 imod1 = 104391568; // modinv(MOD1, MOD2);
    static const mint2 imod01 = 187290749; // imod1 / MOD0;

    // small case (T = mint, long long)
    template<class T> vector<T> naive_mul(const vector<T> &A, const vector<T> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        vector<T> res(N + M - 1);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < M; ++j)
                res[i + j] += A[i] * B[j];
        return res;
    }

    // mul by convolution
    template<class mint> vector<mint> mul(const vector<mint> &A, const vector<mint> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        if (min(N, M) < 30) return naive_mul(A, B);
        int MOD = A[0].get_mod();
        int size_fft = get_fft_size(N, M);
        if (MOD == 998244353) {
            vector<mint> a(size_fft), b(size_fft), c(size_fft);
            for (int i = 0; i < N; ++i) a[i] = A[i];
            for (int i = 0; i < M; ++i) b[i] = B[i];
            trans(a), trans(b);
            vector<mint> res(size_fft);
            for (int i = 0; i < size_fft; ++i) res[i] = a[i] * b[i];
            trans(res, true);
            res.resize(N + M - 1);
            return res;
        }
        vector<mint0> a0(size_fft, 0), b0(size_fft, 0), c0(size_fft, 0);
        vector<mint1> a1(size_fft, 0), b1(size_fft, 0), c1(size_fft, 0);
        vector<mint2> a2(size_fft, 0), b2(size_fft, 0), c2(size_fft, 0);
        for (int i = 0; i < N; ++i)
            a0[i] = A[i].val, a1[i] = A[i].val, a2[i] = A[i].val;
        for (int i = 0; i < M; ++i)
            b0[i] = B[i].val, b1[i] = B[i].val, b2[i] = B[i].val;
        trans(a0), trans(a1), trans(a2), trans(b0), trans(b1), trans(b2);
        for (int i = 0; i < size_fft; ++i) {
            c0[i] = a0[i] * b0[i];
            c1[i] = a1[i] * b1[i];
            c2[i] = a2[i] * b2[i];
        }
        trans(c0, true), trans(c1, true), trans(c2, true);
        mint mod0 = MOD0, mod01 = mod0 * MOD1;
        vector<mint> res(N + M - 1);
        for (int i = 0; i < N + M - 1; ++i) {
            int y0 = c0[i].val;
            int y1 = (imod0 * (c1[i] - y0)).val;
            int y2 = (imod01 * (c2[i] - y0) - imod1 * y1).val;
            res[i] = mod01 * y2 + mod0 * y1 + y0;
        }
        return res;
    }
};

// Formal Power Series
template<typename mint> struct FPS : vector<mint> {
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
        if (this->empty()) this->resize(1);
        (*this)[0] += v;
        return *this;
    }
    constexpr FPS& operator += (const FPS &r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return this->normalize();
    }
    constexpr FPS& operator -= (const mint &v) {
        if (this->empty()) this->resize(1);
        (*this)[0] -= v;
        return *this;
    }
    constexpr FPS& operator -= (const FPS &r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return this->normalize();
    }
    constexpr FPS& operator *= (const mint &v) {
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= v;
        return *this;
    }
    constexpr FPS& operator *= (const FPS &r) {
        return *this = NTT::mul((*this), r);
    }
    constexpr FPS& operator /= (const mint &v) {
        assert(v != 0);
        mint iv = modinv(v);
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= iv;
        return *this;
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
    constexpr mint eval(const mint &v) {
        mint res = 0;
        for (int i = (int)this->size()-1; i >= 0; --i) {
            res *= v;
            res += (*this)[i];
        }
        return res;
    }

    // advanced operation
    // df/dx
    constexpr FPS diff() const {
        int n = (int)this->size();
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
    constexpr FPS inv(int deg) const {
        assert((*this)[0] != 0);
        if (deg < 0) deg = (int)this->size();
        FPS res({mint(1) / (*this)[0]});
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + res - res * res * pre(i << 1)).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS inv() const {
        return inv((int)this->size());
    }
    
    // log(f) = \int f'/f dx, f[0] must be 1
    constexpr FPS log(int deg) const {
        assert((*this)[0] == 1);
        FPS res = (diff() * inv(deg)).integral();
        res.resize(deg);
        return res;
    }
    constexpr FPS log() const {
        return log((int)this->size());
    }
    
    // exp(f), f[0] must be 0
    constexpr FPS exp(int deg) const {
        assert((*this)[0] == 0);
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = res * (pre(i << 1) - res.log(i << 1) + 1).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS exp() const {
        return exp((int)this->size());
    }
    
    // pow(f) = exp(e * log f)
    constexpr FPS pow(long long e, int deg) const {
        if (e == 0) {
            FPS res(deg, 0);
            res[0] = 1;
            return res;
        }
        long long i = 0;
        while (i < (int)this->size() && (*this)[i] == 0) ++i;
        if (i == (int)this->size() || i > (deg - 1) / e) return FPS(deg, 0);
        mint k = (*this)[i];
        FPS res = ((((*this) >> i) / k).log(deg) * e).exp(deg) * mint(k).pow(e) << (e * i);
        res.resize(deg);
        return res;
    }
    constexpr FPS pow(long long e) const {
        return pow(e, (int)this->size());
    }
    
    // sqrt(f), f[0] must be 1
    constexpr FPS sqrt_base(int deg) const {
        assert((*this)[0] == 1);
        mint inv2 = mint(1) / 2;
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + pre(i << 1) * res.inv(i << 1)).pre(i << 1);
            for (mint &x : res) x *= inv2;
        }
        res.resize(deg);
        return res;
    }
    constexpr FPS sqrt_base() const {
        return sqrt_base((int)this->size());
    }
    
    // friend operators
    friend constexpr FPS diff(const FPS &f) { return f.diff(); }
    friend constexpr FPS integral(const FPS &f) { return f.integral(); }
    friend constexpr FPS inv(const FPS &f, int deg) { return f.inv(deg); }
    friend constexpr FPS inv(const FPS &f) { return f.inv((int)f.size()); }
    friend constexpr FPS log(const FPS &f, int deg) { return f.log(deg); }
    friend constexpr FPS log(const FPS &f) { return f.log((int)f.size()); }
    friend constexpr FPS exp(const FPS &f, int deg) { return f.exp(deg); }
    friend constexpr FPS exp(const FPS &f) { return f.exp((int)f.size()); }
    friend constexpr FPS pow(const FPS &f, long long e, int deg) { return f.pow(e, deg); }
    friend constexpr FPS pow(const FPS &f, long long e) { return f.pow(e, (int)f.size()); }
    friend constexpr FPS sqrt_base(const FPS &f, int deg) { return f.sqrt_base(deg); }
    friend constexpr FPS sqrt_base(const FPS &f) { return f.sqrt_base((int)f.size()); }
};


/*/////////////////////////////*/
// Matrix
/*/////////////////////////////*/

// matrix
template<class mint> struct MintMatrix {
    // inner value
    vector<vector<mint>> val;
    
    // constructors
    MintMatrix(int H, int W) : val(H, vector<mint>(W)) {}
    MintMatrix(int H, int W, mint x) : val(H, vector<mint>(W, x)) {}
    MintMatrix(const MintMatrix &mat) : val(mat.val) {}
    void init(int H, int W, mint x) {
        val.assign(H, vector<mint>(W, x));
    }
    void resize(int H, int W) {
        val.resize(H);
        for (int i = 0; i < H; ++i) val[i].resize(W);
    }
    
    // getter and debugger
    constexpr int height() const { return (int)val.size(); }
    constexpr int width() const { return (int)val[0].size(); }
    vector<mint>& operator [] (int i) { return val[i]; }
    constexpr vector<mint>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const MintMatrix<mint> &mat) {
        os << endl;
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) {
                if (j) os << ", ";
                os << mat.val[i][j];
            }
            os << endl;
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
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                val[i][j] = val[i][j] + r.val[i][j];
            }
        }
        return *this;
    }
    constexpr MintMatrix& operator -= (const MintMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                val[i][j] = val[i][j] - r.val[i][j];
            }
        }
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
    constexpr MintMatrix operator + () const { return MintMatrix(*this); }
    constexpr MintMatrix operator - () const { return MintMatrix(*this) *= mint(-1); }
    constexpr MintMatrix operator + (const MintMatrix &r) const { return MintMatrix(*this) += r; }
    constexpr MintMatrix operator - (const MintMatrix &r) const { return MintMatrix(*this) -= r; }
    constexpr MintMatrix operator * (const mint &v) const { return MintMatrix(*this) *= v; }
    constexpr MintMatrix operator * (const MintMatrix &r) const { return MintMatrix(*this) *= r; }
    constexpr vector<mint> operator * (const vector<mint> &v) const {
        assert(width() == v.size());
        vector<mint> res(height(), mint(0));
        for (int i = 0; i < height(); i++)
            for (int j = 0; j < width(); j++)
                res[i] += val[i][j] * v[j];
        return res;
    }
    
    // pow
    constexpr MintMatrix pow(long long n) const {
        assert(height() == width());
        MintMatrix<mint> res(height(), width()), mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = mint(1);
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr MintMatrix<mint> pow(const MintMatrix<mint> &mat, long long n) {
        return mat.pow(n);
    }
    
    // gauss-jordan
    constexpr int find_pivot(int cur_rank, int col) const {
        int pivot = -1;
        for (int row = cur_rank; row < height(); ++row) {
            if (val[row][col] != 0) {
                pivot = row;
                break;
            }
        }
        return pivot;
    }
    constexpr void sweep(int cur_rank, int col, int pivot) {
        swap(val[pivot], val[cur_rank]);
        auto ifac = val[cur_rank][col].inv();
        for (int col2 = 0; col2 < width(); ++col2) {
            val[cur_rank][col2] *= ifac;
        }
        for (int row = 0; row < height(); ++row) {
            if (row != cur_rank && val[row][col] != 0) {
                auto fac = val[row][col];
                for (int col2 = 0; col2 < width(); ++col2) {
                    val[row][col2] -= val[cur_rank][col2] * fac;
                }
            }
        }
    }
    constexpr int gauss_jordan(int not_sweep_width = 0) {
        int rank = 0;
        for (int col = 0; col < width(); ++col) {
            if (col == width() - not_sweep_width) break;
            int pivot = find_pivot(rank, col);
            if (pivot == -1) continue;
            sweep(rank++, col, pivot);
        }
        return rank;
    }
    constexpr int gauss_jordan(int not_sweep_width, vector<int> &core) {
        core.clear();
        int rank = 0;
        for (int col = 0; col < width(); ++col) {
            if (col == width() - not_sweep_width) break;
            int pivot = find_pivot(rank, col);
            if (pivot == -1) continue;
            core.push_back(col);
            sweep(rank++, col, pivot);
        }
        return rank;
    }
    friend constexpr int gauss_jordan(MintMatrix<mint> &mat, int not_sweep_width = 0) {
        return mat.gauss_jordan(not_sweep_width);
    }

    // find one solution
    friend constexpr int linear_equation
    (const MintMatrix<mint> &mat, const vector<mint> &b, vector<mint> &res) {
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
    friend constexpr int linear_equation(const MintMatrix<mint> &mat, const vector<mint> &b) {
        vector<mint> res;
        return linear_equation(mat, b, res);
    }

    // find all solutions
    friend int linear_equation
    (const MintMatrix<mint> &mat, const vector<mint> &b, vector<mint> &res, vector<vector<mint>> &zeros) {
        // extend
        MintMatrix<mint> A(mat.height(), mat.width() + 1);
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) A[i][j] = mat.val[i][j];
            A[i].back() = b[i];
        }
        vector<int> core;
        int rank = A.gauss_jordan(1, core);
        
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
        MintMatrix<mint> A(*this);
        int rank = 0;
        mint res = 1;
        for (int col = 0; col < width(); ++col) {
            int pivot = A.find_pivot(rank, col);
            if (pivot == -1) return mint(0);
            res *= A[pivot][rank];
            A.sweep(rank++, col, pivot);
        }
        return res;
    }
    friend constexpr mint det(const MintMatrix<mint> &mat) {
        return mat.det();
    }
};


/*/////////////////////////////*/
// Flow
/*/////////////////////////////*/

// edge class (for network-flow)
template<class FLOWTYPE> struct FlowEdge {
    // core members
    int rev, from, to;
    FLOWTYPE cap, icap, flow;
    
    // constructor
    FlowEdge(int r, int f, int t, FLOWTYPE c)
    : rev(r), from(f), to(t), cap(c), icap(c), flow(0) {}
    void reset() { cap = icap, flow = 0; }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowEdge& E) {
        return s << E.from << "->" << E.to << '(' << E.flow << '/' << E.icap << ')';
    }
};

// graph class (for network-flow)
template<class FLOWTYPE> struct FlowGraph {
    // core members
    vector<vector<FlowEdge<FLOWTYPE>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    
    // constructor
    FlowGraph(int n = 0) : list(n) { }
    void init(int n = 0) {
        list.assign(n, FlowEdge<FLOWTYPE>());
        pos.clear();
    }
    
    // getter
    vector<FlowEdge<FLOWTYPE>> &operator [] (int i) {
        return list[i];
    }
    const vector<FlowEdge<FLOWTYPE>> &operator [] (int i) const {
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowEdge<FLOWTYPE> &get_rev_edge(const FlowEdge<FLOWTYPE> &e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    FlowEdge<FLOWTYPE> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowEdge<FLOWTYPE> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowEdge<FLOWTYPE>> get_edges() const {
        vector<FlowEdge<FLOWTYPE>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowEdge<FLOWTYPE> &e : list[i]) e.reset();
        }
    }
    void change_edge(FlowEdge<FLOWTYPE> &e, FLOWTYPE new_cap, FLOWTYPE new_flow) {
        FlowEdge<FLOWTYPE> &re = get_rev_edge(e);
        e.cap = new_cap - new_flow, e.icap = new_cap, e.flow = new_flow;
        re.cap = new_flow;
    }
    
    // add_edge
    void add_edge(int from, int to, FLOWTYPE cap) {
        pos.emplace_back(from, (int)list[from].size());
        list[from].push_back(FlowEdge<FLOWTYPE>((int)list[to].size(), from, to, cap));
        list[to].push_back(FlowEdge<FLOWTYPE>((int)list[from].size() - 1, to, from, 0));
    }

    // debug
    friend ostream& operator << (ostream& s, const FlowGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// Dinic
template<class FLOWTYPE> FLOWTYPE Dinic
 (FlowGraph<FLOWTYPE> &G, int s, int t, FLOWTYPE limit_flow)
{
    FLOWTYPE current_flow = 0;
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
            for (const FlowEdge<FLOWTYPE> &e : G[v]) {
                if (level[e.to] < 0 && e.cap > 0) {
                    level[e.to] = level[v] + 1;
                    if (e.to == t) return;
                    que.push(e.to);
                }
            }
        }
    };
    
    // Dinic DFS
    auto dfs = [&](auto self, int v, FLOWTYPE up_flow) {
        if (v == t) return up_flow;
        FLOWTYPE res_flow = 0;
        for (int &i = iter[v]; i < (int)G[v].size(); ++i) {
            FlowEdge<FLOWTYPE> &e = G[v][i], &re = G.get_rev_edge(e);
            if (level[v] >= level[e.to] || e.cap == 0) continue;
            FLOWTYPE flow = self(self, e.to, min(up_flow - res_flow, e.cap));
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
            FLOWTYPE flow = dfs(dfs, s, limit_flow - current_flow);
            if (!flow) break;
            current_flow += flow;
        }
    }
    return current_flow;
};

template<class FLOWTYPE> FLOWTYPE Dinic(FlowGraph<FLOWTYPE> &G, int s, int t) {
    return Dinic(G, s, t, numeric_limits<FLOWTYPE>::max());
}


// edge class (for network-flow)
template<class FLOWTYPE, class COSTTYPE> struct FlowCostEdge {
    // core members
    int rev, from, to;
    FLOWTYPE cap, icap, flow;
    COSTTYPE cost;
    
    // constructor
    FlowCostEdge(int rev, int from, int to, FLOWTYPE cap, COSTTYPE cost)
    : rev(rev), from(from), to(to), cap(cap), icap(cap), flow(0), cost(cost) {}
    void reset() { cap = icap, flow = 0; }
    
    // debug
    friend ostream& operator << (ostream& s, const FlowCostEdge& e) {
        return s << e.from << " -> " << e.to
        << " (" << e.flow << "/" << e.icap << ", " << e.cost << ")";
    }
};

// graph class (for network-flow)
template<class FLOWTYPE, class COSTTYPE> struct FlowCostGraph {
    // core members
    vector<vector<FlowCostEdge<FLOWTYPE, COSTTYPE>>> list;
    vector<pair<int,int>> pos;  // pos[i] := {vertex, order of list[vertex]} of i-th edge
    
    // constructor
    FlowCostGraph(int n = 0) : list(n) { }
    void init(int n = 0) {
        list.assign(n, FlowCostEdge<FLOWTYPE, COSTTYPE>());
        pos.clear();
    }
    
    // getter
    vector<FlowCostEdge<FLOWTYPE, COSTTYPE>> &operator [] (int i) {
        return list[i];
    }
    const vector<FlowCostEdge<FLOWTYPE, COSTTYPE>> &operator [] (int i) const {
        return list[i];
    }
    size_t size() const {
        return list.size();
    }
    FlowCostEdge<FLOWTYPE, COSTTYPE> &get_rev_edge
    (const FlowCostEdge<FLOWTYPE, COSTTYPE> &e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
    FlowCostEdge<FLOWTYPE, COSTTYPE> &get_edge(int i) {
        return list[pos[i].first][pos[i].second];
    }
    const FlowCostEdge<FLOWTYPE, COSTTYPE> &get_edge(int i) const {
        return list[pos[i].first][pos[i].second];
    }
    vector<FlowCostEdge<FLOWTYPE, COSTTYPE>> get_edges() const {
        vector<FlowCostEdge<FLOWTYPE, COSTTYPE>> edges;
        for (int i = 0; i < (int)pos.size(); ++i) {
            edges.push_back(get_edge(i));
        }
        return edges;
    }
    
    // change edges
    void reset() {
        for (int i = 0; i < (int)list.size(); ++i) {
            for (FlowCostEdge<FLOWTYPE, COSTTYPE> &e : list[i]) e.reset();
        }
    }
    
    // add_edge
    void add_edge(int from, int to, FLOWTYPE cap, COSTTYPE cost) {
        pos.emplace_back(from, (int)list[from].size());
        list[from].push_back(FlowCostEdge<FLOWTYPE, COSTTYPE>
                             ((int)list[to].size(), from, to, cap, cost));
        list[to].push_back(FlowCostEdge<FLOWTYPE, COSTTYPE>
                           ((int)list[from].size() - 1, to, from, 0, -cost));
    }

    // debug
    friend ostream& operator << (ostream& s, const FlowCostGraph &G) {
        const auto &edges = G.get_edges();
        for (const auto &e : edges) s << e << endl;
        return s;
    }
};

// min-cost max-flow (<= limit_flow), slope ver.
template<class FLOWTYPE, class COSTTYPE> vector<pair<FLOWTYPE, COSTTYPE>>
MinCostFlowSlope(FlowCostGraph<FLOWTYPE, COSTTYPE> &G, int S, int T, FLOWTYPE limit_flow)
{
    // result values
    FLOWTYPE cur_flow = 0;
    COSTTYPE cur_cost = 0, pre_cost = -1;
    vector<pair<FLOWTYPE, COSTTYPE>> res;
    res.emplace_back(cur_flow, cur_cost);
    
    // intermediate values
    vector<bool> seen((int)G.size(), false);
    vector<COSTTYPE> dist((int)G.size(), numeric_limits<COSTTYPE>::max());
    vector<int> prevv((int)G.size(), -1), preve((int)G.size(), -1);
    
    // dual
    auto dual_step = [&]() -> bool {
        seen.assign((int)G.size(), false);
        dist.assign((int)G.size(), numeric_limits<COSTTYPE>::max());
        seen[S] = true, dist[S] = 0;
        while (true) {
            bool update = false;
            for (int v = 0; v < (int)G.size(); ++v) {
                if (!seen[v]) continue;
                for (int i = 0; i < G[v].size(); ++i) {
                    const FlowCostEdge<FLOWTYPE, COSTTYPE> &e = G[v][i];
                    if (e.cap > 0 && (!seen[e.to] || dist[e.to] > dist[v] + e.cost)) {
                        dist[e.to] = dist[v] + e.cost;
                        prevv[e.to] = v;
                        preve[e.to] = i;
                        seen[e.to] = true;
                        update = true;
                    }
                }
            }
            if (!update) break;
        }
        return seen[T];
    };
    
    // primal
    auto primal_step = [&]() -> void {
        FLOWTYPE flow = limit_flow - cur_flow;
        COSTTYPE cost = dist[T];
        for (int v = T; v != S; v = prevv[v]) {
            flow = min(flow, G[prevv[v]][preve[v]].cap);
        }
        for (int v = T; v != S; v = prevv[v]) {
            FlowCostEdge<FLOWTYPE, COSTTYPE> &e = G[prevv[v]][preve[v]];
            FlowCostEdge<FLOWTYPE, COSTTYPE> &re = G.get_rev_edge(e);
            e.cap -= flow, e.flow += flow;
            re.cap += flow, re.flow -= flow;
        }
        cur_flow += flow;
        cur_cost += flow * cost;
        if (pre_cost == cost) res.pop_back();
        res.emplace_back(cur_flow, cur_cost);
        pre_cost = cur_cost;
    };
    
    // primal-dual
    while (cur_flow < limit_flow) {
        if (!dual_step()) break;
        primal_step();
    }
    return res;
}

// min-cost max-flow, slope ver.
template<class FLOWTYPE, class COSTTYPE> vector<pair<FLOWTYPE, COSTTYPE>>
MinCostFlowSlope(FlowCostGraph<FLOWTYPE, COSTTYPE> &G, int S, int T)
{
    return MinCostFlowSlope(G, S, T, numeric_limits<FLOWTYPE>::max());
}

// min-cost max-flow (<= limit_flow)
template<class FLOWTYPE, class COSTTYPE> pair<FLOWTYPE, COSTTYPE>
MinCostFlow(FlowCostGraph<FLOWTYPE, COSTTYPE> &G, int S, int T, FLOWTYPE limit_flow)
{
    return MinCostFlowSlope(G, S, T, limit_flow).back();
}

// min-cost max-flow (<= limit_flow)
template<class FLOWTYPE, class COSTTYPE> pair<FLOWTYPE, COSTTYPE>
MinCostFlow(FlowCostGraph<FLOWTYPE, COSTTYPE> &G, int S, int T)
{
    return MinCostFlow(G, S, T, numeric_limits<FLOWTYPE>::max());
}


// Hopcroft-Karp
struct HopcroftKarp {
    // input
    int size_left, size_right;
    vector<vector<int>> list; // left to right

    // results
    vector<int> lr, rl;
    
    // intermediate results
    vector<bool> seen, matched;
    vector<int> level;
    
    // constructor
    HopcroftKarp(int l, int r) : size_left(l), size_right(r), list(l, vector<int>()) { }
    void add_edge(int from, int to) {
        list[from].push_back(to);
    }

    // getter, debugger
    const vector<int> &operator [] (int i) const { 
        return list[i];
    }
    friend ostream& operator << (ostream& s, const HopcroftKarp& G) {
        s << endl;
        for (int i = 0; i < G.list.size(); ++i) {
            s << i << ": ";
            for (int j = 0; j < G.list[i].size(); ++j) {
                s << G.list[i][j];
                if (j + 1 != G.list[i].size()) s << ", ";
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
            for (int i = 0; i < list[left].size(); ++i) {
                int right = list[left][i];
                int next = rl[right];
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
        for (int i = 0; i < list[left].size(); ++i) {
            int right = list[left][i];
            int next = rl[right];
            if (next == -1) next = size_left;
            if (level[next] > level[left] && hodfs(next)) {
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
            if (rl[r] != -1) lr[rl[r]] = r;
        }
        return res;
    }
};


/*/////////////////////////////*/
// Tree
/*/////////////////////////////*/

// Run Tree (including Euler Tour)
template<class Graph = vector<vector<int>>> struct RunTree {
    // id[v][w] := the index of node w in G[v]
    vector<unordered_map<int, int>> id;

    // num[v][i] := the size of subtree of G[v][i] with parent v
    vector<vector<long long>> num;
    
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
    RunTree(const Graph &G, int root = 0) : root(root) {
        init(G, root);
    }
    
    // init
    void init(const Graph &G, int root = 0) {
        int N = (int)G.size();
        id.assign(N, unordered_map<int,int>()), num.assign(N, vector<long long>());
        for (int v = 0; v < N; v++) num[v].assign((int)G[v].size(), 0);
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
        return num[u][id[u][v]];
    }

    // get first / last id of node v in Euler tour
    int vs(int v) { return v_s_id[v]; }
    int vt(int v) { return v_t_id[v]; }

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

    // get_parent(v, p) := the parent of v directed for p
    int get_parent(int v, int p) {
        if (v == p) return -1;
        int lca = get_lca(v, p);
        if (lca != v) return parent[0][v];
        for (int i = (int)parent.size()-1; i >= 0; i--) {
            if (parent[i][p] != -1 && depth[parent[i][p]] > depth[v]) {
                p = parent[i][p];
            }
        }
        return p;
    }
    
    // rec
    int rec(const Graph &G, int v, int p, int d, int &ord) {
        int p_index = -1;
        int sum = 1;
        parent[0][v] = p, depth[v] = d;
        tour[ord] = v, v_s_id[v] = v_t_id[v] = ord;
        ord++;
        for (int i = 0; i < (int)G[v].size(); i++) {
            int ch = G[v][i];
            id[v][ch] = i;
            if (ch == p) {
                p_index = i;
                continue;
            }
            e_id[ch * 2] = ord - 1;
            int s = rec(G, ch, v, d+1, ord);
            num[v][i] = s;
            sum += s;
            tour[ord] = v;
            v_t_id[v] = ord;
            e_id[ch * 2 + 1] = ord - 1;
            ord++;
        }
        if (p_index != -1) num[v][p_index] = (int)G.size() - sum;
        return sum;
    }
};


/*/////////////////////////////*/
// Union-Find
/*/////////////////////////////*/

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
    
    bool merge(int x, int y, bool merge_technique = true) {
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


/*/////////////////////////////////*/
// Sparse Table, Binary Indexed Tree
/*/////////////////////////////////*/

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
            s << x << " ";
        }
        return s;
    }
};

// Sparse Table
template<class MeetSemiLattice> struct SparseTable {
    using Func = function<MeetSemiLattice(MeetSemiLattice, MeetSemiLattice)>;

    // core member
    Func OP;
    vector<vector<MeetSemiLattice>> dat;
    vector<int> height;
    
    SparseTable() { }
    SparseTable(const vector<MeetSemiLattice> &vec, const Func &op)  { init(vec, op); }
    void init(const vector<MeetSemiLattice> &vec, const Func &op) {
        OP = op;
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
    
    MeetSemiLattice get(int a, int b) {
        return OP(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};


/*/////////////////////////////*/
// Segment Tree
/*/////////////////////////////*/

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
template<class T> auto seg_op_min = [](T p, T q) { return min(p, q); };
template<class T> auto seg_op_max = [](T p, T q) { return max(p, q); };
template<class T> auto seg_op_min_with_index = [](pair<T,int> p, pair<T,int> q) { return min(p, q); };
template<class T> auto seg_op_max_with_index = [](pair<T,int> p, pair<T,int> q) { return max(p, q); };
template<class T> SegmentTree<T> RangeMin() {
    return SegmentTree<T>(seg_op_min<T>, numeric_limits<T>::max()/2);
}
template<class Monoid> SegmentTree<Monoid> RangeMax() {
    return SegmentTree<Monoid>(seg_op_max<Monoid>, -numeric_limits<Monoid>::max()/2);
}
template<class Monoid> SegmentTree<pair<Monoid,int>> RangeMinWithIndex() {
    return SegmentTree<Monoid>(seg_op_min_with_index<Monoid>, {numeric_limits<Monoid>::max()/2, -1});
}
template<class Monoid> SegmentTree<pair<Monoid,int>> RangeMaxWithIndex() {
    return SegmentTree<Monoid>(seg_op_max_with_index<Monoid>, {-numeric_limits<Monoid>::max()/2, -1});
}

// various lazy segment trees
template<class Monoid> auto seg_op_add = [](Monoid p, Monoid q) { return p + q; };
template<class Monoid> auto seg_op_add_with_sum = [](pair<Monoid,long long> p, pair<Monoid,long long> q) { 
    return make_pair(p.first + q.first, p.second + q.second);
};
template<class Monoid, class Action> auto seg_act_change = [](Action f, Monoid x) {
    return (f != -1 ? f : x);
};
template<class Monoid, class Action> auto seg_act_add = [](Action f, Monoid x) {
    return f + x;
};
template<class Monoid, class Action> auto seg_act_change_with_index = [](Action f, pair<Monoid,int> x) {
    return (f != -1 ? make_pair(f, x.second) : x);
};
template<class Monoid, class Action> auto seg_act_add_with_index = [](Action f, pair<Monoid,int> x) {
    return make_pair(f.first + x, f.second);
};
template<class Monoid, class Action> auto seg_act_change_with_sum = [](Action f, Monoid x) {
    return (f != -1 ? make_pair(f * x.second, x.second) : x);
};
template<class Action> auto seg_comp_change = [](Action g, Action f) {
    return (g != -1 ? g : f);
};
template<class Action> auto seg_comp_add = [](Action g, Action f) {
    return g + f;
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMin() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_min<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        numeric_limits<Monoid>::max()/2, -1);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMax() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_max<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        -numeric_limits<Monoid>::max()/2, -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,int>, Action> RangeChangeRangeMinWithIndex() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_min_with_index<Monoid>, seg_act_change_with_index<Monoid, Action>, seg_comp_change<Action>,
        make_pair(numeric_limits<Monoid>::max()/2 -1), -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,int>, Action> RangeChangeRangeMaxWithIndex() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_max_with_index<Monoid>, seg_act_change_with_index<Monoid, Action>, seg_comp_change<Action>,
        make_pair(-numeric_limits<Monoid>::max()/2 -1), -1);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeAddRangeMin() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_min<Monoid>, seg_act_add<Monoid, Action>, seg_comp_add<Action>,
        numeric_limits<Monoid>::max()/2, -1);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeAddRangeMax() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_max<Monoid>, seg_act_add<Monoid, Action>, seg_comp_add<Action>,
        -numeric_limits<Monoid>::max()/2, -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,int>, Action> RangeAddRangeMinWithIndex() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_min_with_index<Monoid>, seg_act_add_with_index<Monoid, Action>, seg_comp_add<Action>,
        make_pair(numeric_limits<Monoid>::max()/2, -1), -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,int>, Action> RangeAddRangeMaxWithIndex() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_max_with_index<Monoid>, seg_act_add_with_index<Monoid, Action>, seg_comp_add<Action>,
        make_pair(-numeric_limits<Monoid>::max()/2, -1), -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,long long>, Action> RangeChangeRangeSum() {
    return LazySegmentTree<Monoid, Action>(
        seg_op_add_with_sum<Monoid>, seg_act_change_with_sum<Monoid, Action>, seg_comp_change<Action>,
        make_pair(Monoid(0), 1), -1);
};



/*/////////////////////////////*/
// Solver
/*/////////////////////////////*/

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    // FastRead Read
    // FastWrite Write
    
    
}