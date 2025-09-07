// code template is in https://github.com/drken1215/algorithm/blob/master/template_atcoder.cpp
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Utility
//------------------------------//

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

// 4-neighbor
const vector<int> dx = {1, 0, -1, 0};
const vector<int> dy = {0, 1, 0, -1};

// 8-neighbor
const vector<int> dx8 = {1, 0, -1, 0, 1, -1, 1, -1};
const vector<int> dy8 = {0, 1, 0, -1, 1, 1, -1, -1};

// min non-negative i such that n <= 2^i
int ceil_pow2(int n) {
    int i = 0;
    while ((1U << i) < (unsigned int)(n)) i++;
    return i;
}

// num of i such that (x & (1 << i)) != 0
int popcnt(int x) { return __builtin_popcount(x); }
int popcnt(unsigned int x) { return __builtin_popcount(x); }
int popcnt(long long x) { return __builtin_popcountll(x); }
int popcnt(unsigned long long x) { return __builtin_popcountll(x); }

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
long long randInt(long long minv, long long maxv) {
    long long a = randInt(), b = randInt();
    return (a * (1LL<<29) + b) % (maxv - minv + 1) + minv;
}
template<class T> void shuffle(vector<T>& vec) {
    int n = vec.size();
    for (int i = n - 1; i > 0; --i) {
        int k = randInt() % (i + 1);
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
u128 to_uinteger(const string &s) {
    u128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    return res;
}
istream& operator >> (istream &is, u128 &x) {
    string s;
    is >> s;
    x = to_uinteger(s);
    return is;
}
ostream& operator << (ostream &os, const u128 &x) {
    u128 ax = x;
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


//------------------------------//
// Fast IO
//------------------------------//

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
    template<class T> static constexpr int DIGITS = numeric_limits<T>::digits10 + 1;
    template<class T> static constexpr auto POW10 = [] {
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


//------------------------------//
// mod algorithms
//------------------------------//

// safe mod
template<class T_VAL, class T_MOD>
constexpr T_VAL safe_mod(T_VAL a, T_MOD m) {
    assert(m > 0);
    a %= m;
    if (a < 0) a += m;
    return a;
}

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

// mod inv
template<class T_VAL, class T_MOD>
constexpr T_VAL mod_inv(T_VAL a, T_MOD m) {
    T_VAL b = m, u = 1, v = 0;
    while (b > 0) {
        T_VAL t = a / b;
        a -= t * b, swap(a, b);
        u -= t * v, swap(u, v);
    }
    u %= m;
    if (u < 0) u += m;
    return u;
}

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
        if (PRIME) {
            assert(val);
            return pow(get_umod() - 2);
        } else {
            assert(val);
            return mod_inv((long long)(val), get_umod());
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
    friend constexpr istream& operator >> (istream &is, Fp<MOD> &x) {
        long long tmp = 1;
        is >> tmp;
        tmp = tmp % (long long)(get_umod());
        if (tmp < 0) tmp += get_umod();
        x.val = (unsigned int)(tmp);
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
        return mod_inv((long long)(val), get_umod());
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


// Garner's algorithm
// for each step, we solve "coeffs[k] * t[k] + constants[k] = b[k] (mod. m[k])"
//      coeffs[k] = m[0]m[1]...m[k-1]
//      constants[k] = t[0] + t[1]m[0] + ... + t[k-1]m[0]m[1]...m[k-2]

// if m is not coprime, call this function first
template<class T_VAL>
bool preGarner(vector<T_VAL> &b, vector<T_VAL> &m) {
    assert(b.size() == m.size());
    T_VAL res = 1;
    for (int i = 0; i < (int)b.size(); i++) {
        for (int j = 0; j < i; ++j) {
            T_VAL g = gcd(m[i], m[j]);
            if ((b[i] - b[j]) % g != 0) return false;
            m[i] /= g, m[j] /= g;
            T_VAL gi = gcd(m[i], g), gj = g/gi;
            do {
                g = gcd(gi, gj);
                gi *= g, gj /= g;
            } while (g != 1);
            m[i] *= gi, m[j] *= gj;
            b[i] %= m[i], b[j] %= m[j];
        }
    }
    vector<T_VAL> b2, m2;
    for (int i = 0; i < (int)b.size(); i++) {
        if (m[i] == 1) continue;
        b2.emplace_back(b[i]), m2.emplace_back(m[i]);
    }
    b = b2, m = m2;
    return true;
}

// find x (%MOD), LCM (%MOD) (m must be coprime)
template<class T_VAL>
T_VAL Garner(vector<T_VAL> b, vector<T_VAL> m) {
    assert(b.size() == m.size());
    using mint = DynamicModint;
    int num = (int)m.size();
    T_VAL res = 0, lcm = 1;
    vector<long long> coeffs(num, 1), constants(num, 0);
    for (int k = 0; k < num; k++) {
        mint::set_mod(m[k]);
        T_VAL t = ((mint(b[k]) - constants[k]) / coeffs[k]).val;
        for (int i = k + 1; i < num; i++) {
            constants[i] = safe_mod(constants[i] + t * coeffs[i], m[i]);
            coeffs[i] = safe_mod(coeffs[i] * m[k], m[i]);
        }
        res += t * lcm;
        lcm *= m[k];
    }
    return res;
}

// find x, LCM (m must be coprime)
template<class T_VAL, class T_MOD>
T_VAL Garner(vector<T_VAL> b, vector<T_VAL> m, T_MOD MOD) {
    assert(b.size() == m.size());
    assert(MOD > 0);
    using mint = DynamicModint;
    int num = (int)m.size();
    T_VAL res = 0, lcm = 1;
    vector<long long> coeffs(num, 1), constants(num, 0);
    for (int k = 0; k < num; k++) {
        mint::set_mod(m[k]);
        T_VAL t = ((mint(b[k]) - constants[k]) / coeffs[k]).val;
        for (int i = k + 1; i < num; i++) {
            constants[i] = safe_mod(constants[i] + t * coeffs[i], m[i]);
            coeffs[i] = safe_mod(coeffs[i] * m[k], m[i]);
        }
        res = safe_mod(res + t * lcm, MOD);
        lcm = safe_mod(lcm * m[k], MOD);
    }
    return res;
}


//------------------------------//
// Prime
//------------------------------//

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

// montgomery modint (MOD < 2^62, MOD is odd)
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
    mint inv() const { return pow(MOD - 2); }
    mint pow(u128 n) const {
        mint res(1), mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
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

// various methods mod prime P
struct PrimeProcessor {
    using mint = MontgomeryModInt64;
    
    // input prime
    long long prime;
    vector<pair<long long, long long>> pf;  // prime factorization of p-1
    
    // constructors
    PrimeProcessor() {}
    PrimeProcessor(long long p) : prime(p) {
        init(p);
    }
    
    // initializer
    void init(long long p) {
        assert(is_prime(p));
        prime = p;
        if (p % 2 == 1) {
            assert(p < (1LL<<62));
            prime = p;
            pf = prime_factorize(prime - 1);
            mint::set_mod(prime);
        }
    }
    
    // min: x s.t. a^x \equiv 1 (mod prime)
    long long calc_order(long long a) {
        assert(a != 0);
        if (prime == 2) return 1;
        long long res = prime - 1;
        for (const auto &[p, num] : pf) {
            while (res % p == 0 && mint(a).pow(res / p) == 1) res /= p;
        }
        return res;
    }
};


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
    int h = ceil_pow2(n);
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
                    rot *= setup.rate2[bsf(~(unsigned int)(s))];
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
                    rot *= setup.rate3[bsf(~(unsigned int)(s))];
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
    int h = ceil_pow2(n);
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
                    irot *= setup.irate2[bsf(~(unsigned int)(s))];
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
                    irot *= setup.irate3[bsf(~(unsigned int)(s))];
                }
            }
            len -= 2;
        }
    }
    mint in = mint(n).inv();
    for (int i = 0; i < n; i++) v[i] *= in;
}

// naive convolution
template<class T>
vector<T> sub_convolution_naive(const vector<T> &a, const vector<T> &b) {
    int n = (int)a.size(), m = (int)b.size();
    vector<T> res(n + m - 1);
    if (n < m) {
        for (int j = 0; j < m; j++) for (int i = 0; i < n; i++) res[i + j] += a[i] * b[j];
    } else {
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) res[i + j] += a[i] * b[j];
    }
    return res;
}

// ntt convolution
template<class mint>
vector<mint> sub_convolution_ntt(vector<mint> a, vector<mint> b) {
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

// convolution in general mod
template<class mint>
vector<mint> convolution(const vector<mint> &a, const vector<mint> &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    if (min(n, m) <= 60) return sub_convolution_naive(std::move(a), std::move(b));
    if constexpr (std::is_same_v<mint, Fp<998244353>>) return sub_convolution_ntt(a, b);

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
    auto c0 = sub_convolution_ntt(std::move(a0), std::move(b0));
    auto c1 = sub_convolution_ntt(std::move(a1), std::move(b1));
    auto c2 = sub_convolution_ntt(std::move(a2), std::move(b2));

    vector<mint> res(n + m - 1);
    mint mod0 = MOD0, mod01 = mod0 * MOD1;
    for (int i = 0; i < n + m - 1; ++i) {
        unsigned int y0 = c0[i].val;
        unsigned int y1 = (imod0 * (c1[i] - y0)).val;
        unsigned int y2 = (imod01 * (c2[i] - y0) - imod1 * y1).val;
        res[i] = mod01 * y2 + mod0 * y1 + y0;
    }
    return res;
}

// convolution long long (especially, mod 2^64)
vector<unsigned long long> convolution_ull(const vector<unsigned long long> &a, const vector<unsigned long long> &b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    if (min(n, m) <= 60) return sub_convolution_naive(std::move(a), std::move(b));

    static constexpr int MOD0 = 754974721;  // 2^24
    static constexpr int MOD1 = 167772161;  // 2^25
    static constexpr int MOD2 = 469762049;  // 2^26
    static constexpr int MOD3 = 998244353;  // 2^23
    static constexpr int MOD4 = 645922817;  // 2^23
    static constexpr int MOD5 = 897581057;  // 2^23
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    using mint3 = Fp<MOD3>;
    using mint4 = Fp<MOD4>;
    using mint5 = Fp<MOD5>;

    vector<mint0> a0(n, 0), b0(m, 0);
    vector<mint1> a1(n, 0), b1(m, 0);
    vector<mint2> a2(n, 0), b2(m, 0);
    vector<mint3> a3(n, 0), b3(m, 0);
    vector<mint4> a4(n, 0), b4(m, 0);
    vector<mint5> a5(n, 0), b5(m, 0);
    for (int i = 0; i < n; ++i) {
        a0[i] = a[i] % MOD0;
        a1[i] = a[i] % MOD1;
        a2[i] = a[i] % MOD2;
        a3[i] = a[i] % MOD3;
        a4[i] = a[i] % MOD4;
        a5[i] = a[i] % MOD5;
    }
    for (int i = 0; i < m; ++i) {
        b0[i] = b[i] % MOD0;
        b1[i] = b[i] % MOD1;
        b2[i] = b[i] % MOD2;
        b3[i] = b[i] % MOD3;
        b4[i] = b[i] % MOD4;
        b5[i] = b[i] % MOD5;
    }
    auto c0 = sub_convolution_ntt(std::move(a0), std::move(b0));
    auto c1 = sub_convolution_ntt(std::move(a1), std::move(b1));
    auto c2 = sub_convolution_ntt(std::move(a2), std::move(b2));
    auto c3 = sub_convolution_ntt(std::move(a3), std::move(b3));
    auto c4 = sub_convolution_ntt(std::move(a4), std::move(b4));
    auto c5 = sub_convolution_ntt(std::move(a5), std::move(b5));

    vector<unsigned long long> res(n + m - 1);
    for (int i = 0; i < n + m - 1; i++) {
        vector<unsigned long long> rems = {c0[i].val, c1[i].val, c2[i].val, c3[i].val, c4[i].val, c5[i].val};
        vector<unsigned long long> mods = {MOD0, MOD1, MOD2, MOD3, MOD4, MOD5};
        res[i] = Garner(rems, mods);
    }
    return res;
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
    constexpr FPS& operator /= (const mint &v) {
        assert(v != 0);
        mint iv = v.inv();
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

    // polynomial taylor shift
    constexpr FPS taylor_shift(long long c) const {
        int N = (int)this->size() - 1;
        BiCoef<mint> bc(N + 1);
        FPS<mint> p(N + 1), q(N + 1);
        for (int i = 0; i <= N; i++) {
            p[i] = (*this)[i] * bc.fact(i);
            q[N - i] = mint(c).pow(i) * bc.finv(i);
        }
        FPS<mint> pq = p * q;
        FPS<mint> res(N + 1);
        for (int i = 0; i <= N; i++) res[i] = pq[i + N] * bc.finv(i);
        return res;
    }
    
    // friend operators
    friend constexpr FPS diff(const FPS &f) { return f.diff(); }
    friend constexpr FPS integral(const FPS &f) { return f.integral(); }
    friend constexpr FPS inv(const FPS &f, int deg = -1) { return f.inv(deg); }
    friend constexpr FPS log(const FPS &f, int deg = -1) { return f.log(deg); }
    friend constexpr FPS exp(const FPS &f, int deg = -1) { return f.exp(deg); }
    friend constexpr FPS pow(const FPS &f, long long e, int deg = -1) { return f.pow(e, deg); }
    friend constexpr FPS sqrt(const FPS &f, int deg = -1) { return f.sqrt(deg); }
    friend constexpr FPS taylor_shift(const FPS &f, long long c) { return f.taylor_shift(c); }
};

// Bostan-Mori
// find [x^N] P(x)/Q(x), O(K log K log N)
// deg(Q(x)) = K, deg(P(x)) < K
template<typename mint> mint BostanMori(const FPS<mint> &P, const FPS<mint> &Q, long long N) {
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
    MintMatrix() : H(0), W(0) {}
    MintMatrix(int h, int w) : H(h), W(w), val(h, vector<mint>(w)) {}
    MintMatrix(int h, int w, mint x) : H(h), W(w), val(h, vector<mint>(w, x)) {}
    MintMatrix(const MintMatrix &mat) : H(mat.H), W(mat.W), val(mat.val) {}
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
    constexpr vector<mint>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const MintMatrix<mint> &mat) {
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

    // transpose
    constexpr MintMatrix trans() const {
        MintMatrix<mint> res(width(), height());
        for (int row = 0; row < width(); row++) for (int col = 0; col < height(); col++) {
            res[row][col] = val[col][row];
        }
        return res;
    }
    friend constexpr MintMatrix<mint> trans(const MintMatrix<mint> &mat) {
        return mat.trans();
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
    friend constexpr int gauss_jordan(MintMatrix<mint> &mat, int not_sweep_width = 0, bool sweep_upper = true) {
        return mat.gauss_jordan(not_sweep_width, sweep_upper);
    }

    // rank
    constexpr int get_rank() const {
        if (height() == 0 || width() == 0) return 0;
        MintMatrix A(*this);
        if (height() < width()) A = A.trans();
        return A.gauss_jordan(0, false);
    }
    friend constexpr int get_rank(const MintMatrix<mint> &mat) {
        return mat.get_rank();
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
    friend constexpr mint det(const MintMatrix<mint> &mat) {
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
    friend constexpr mint det_nonprime_mod(const MintMatrix<mint> &mat) {
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
        for (int i = 0; i < height(); ++i) for (int j = 0; j < width(); ++j) res[i][j] = A[i][j+width()];
        return res;
    }
    friend constexpr MintMatrix<mint> inv(const MintMatrix<mint> &mat) {
        return mat.inv();
    }
};

// characteristic_polynomial
// find f(x) = det(xI - B), for N x N mint matrix B, O(N^3)
template<class mint>
void hessenberg_reduction(MintMatrix<mint> &M) {
    assert(M.height() == M.width());
    int N = (int)M.height();
    for (int r = 0; r < N - 2; r++) {
        int pivot = -1;
        for (int h = r + 1; h < N; h++) {
            if (M[h][r] != 0) {
                pivot = h;
                break;
            }
        }
        if (pivot == -1) continue;
        for (int i = 0; i < N; i++) swap(M[r + 1][i], M[pivot][i]);
        for (int i = 0; i < N; i++) swap(M[i][r + 1], M[i][pivot]);
        mint ir = M[r + 1][r].inv();
        for (int i = r + 2; i < N; i++) {
            mint ir2 = M[i][r] * ir;
            for (int j = 0; j < N; j++) M[i][j] -= M[r + 1][j] * ir2;
            for (int j = 0; j < N; j++) M[j][r + 1] += M[j][i] * ir2;
        }
    }
}
template<class mint>
FPS<mint> calc_characteristic_polynomial(MintMatrix<mint> B) {
    assert(B.height() == B.width());
    hessenberg_reduction(B);
    int N = (int)B.height();
    vector<FPS<mint>> res(N + 1);
    res[0] = {mint(1)};
    for (int i = 0; i < N; i++) {
        res[i + 1].assign(i + 2, mint(0));
        for (int j = 0; j < i + 1; j++) res[i + 1][j + 1] += res[i][j];
        for (int j = 0; j < i + 1; j++) res[i + 1][j] -= res[i][j] * B[i][i];
        mint beta = 1;
        for (int j = i - 1; j >= 0; j--) {
            beta *= B[j + 1][j];
            mint beta2 = -B[j][i] * beta;
            for (int k = 0; k < j + 1; k++) res[i + 1][k] += beta2 * res[j][k];
        }
    }
    return res[N];
}

// find f(x) = det(M1x + M0), for given N x N matrix M0, M1, in O(N^3)
template<class mint>
FPS<mint> calc_det_linear_expression(MintMatrix<mint> M0, MintMatrix<mint> M1) {
    int N = (int)M0.height();
    assert(M0.width() == N && M1.height() == N && M1.width() == N);
    int con = 0;
    mint invAB = 1;
    for (int p = 0; p < N; p++) {
        int pivot = -1;
        for (int r = p; r < N; r++) {
            if (M1[r][p] != mint(0)) {
                pivot = r;
                break;
            }
        }
        if (pivot == -1) {
            con++;
            if (con > N) return FPS<mint>(N + 1, mint(0));
            for (int r = 0; r < p; r++) {
                mint v = M1[r][p];
                M1[r][p] = 0;
                for (int i = 0; i < N; i++) M0[i][p] -= v * M0[i][r];
            }
            for (int i = 0; i < N; i++) swap(M0[i][p], M1[i][p]);
            p--;
            continue;
        }
        if (pivot != p) {
            swap(M1[p], M1[pivot]);
            swap(M0[p], M0[pivot]);
            invAB *= -1;
        }
        mint v = M1[p][p], iv = v.inv();
        invAB *= v;
        for (int c = 0; c < N; c++) M0[p][c] *= iv, M1[p][c] *= iv;
        for (int r = 0; r < N; r++) {
            if (r == p) continue;
            mint v = M1[r][p];
            for (int c = 0; c < N; c++) M0[r][c] -= M0[p][c] * v, M1[r][c] -= M1[p][c] * v;
        }
    }
    for (int r = 0; r < M0.height(); r++) for (int c = 0; c < M0.width(); c++) M0[r][c] *= -1;
    auto pol = calc_characteristic_polynomial(M0);
    for (auto &x : pol) x *= invAB;
    pol.erase(pol.begin(), pol.begin() + con);
    pol.resize(N + 1);
    return pol;
}


//------------------------------//
// Flow
//------------------------------//

// edge class (for max-flow)
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

// graph class (for max-flow)
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
template<class T> SegmentTree<T> RangeMin(int N = 0) {
    return SegmentTree<T>(N, seg_op_min<T>, numeric_limits<T>::max()/2);
}
template<class Monoid> SegmentTree<Monoid> RangeMax(int N = 0) {
    return SegmentTree<Monoid>(N, seg_op_max<Monoid>, -numeric_limits<Monoid>::max()/2);
}
template<class Monoid> SegmentTree<pair<Monoid,int>> RangeMinWithIndex(int N = 0) {
    return SegmentTree<Monoid>(N, seg_op_min_with_index<Monoid>, {numeric_limits<Monoid>::max()/2, -1});
}
template<class Monoid> SegmentTree<pair<Monoid,int>> RangeMaxWithIndex(int N = 0) {
    return SegmentTree<Monoid>(N, seg_op_max_with_index<Monoid>, {-numeric_limits<Monoid>::max()/2, -1});
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
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMin(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_min<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        numeric_limits<Monoid>::max()/2, -1);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeChangeRangeMax(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_max<Monoid>, seg_act_change<Monoid, Action>, seg_comp_change<Action>,
        -numeric_limits<Monoid>::max()/2, -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,int>, Action> RangeChangeRangeMinWithIndex(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_min_with_index<Monoid>, seg_act_change_with_index<Monoid, Action>, seg_comp_change<Action>,
        make_pair(numeric_limits<Monoid>::max()/2 -1), -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,int>, Action> RangeChangeRangeMaxWithIndex(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_max_with_index<Monoid>, seg_act_change_with_index<Monoid, Action>, seg_comp_change<Action>,
        make_pair(-numeric_limits<Monoid>::max()/2 -1), -1);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeAddRangeMin(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_min<Monoid>, seg_act_add<Monoid, Action>, seg_comp_add<Action>,
        numeric_limits<Monoid>::max()/2, -1);
};
template<class Monoid, class Action> LazySegmentTree<Monoid, Action> RangeAddRangeMax(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_max<Monoid>, seg_act_add<Monoid, Action>, seg_comp_add<Action>,
        -numeric_limits<Monoid>::max()/2, -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,int>, Action> RangeAddRangeMinWithIndex(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_min_with_index<Monoid>, seg_act_add_with_index<Monoid, Action>, seg_comp_add<Action>,
        make_pair(numeric_limits<Monoid>::max()/2, -1), -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,int>, Action> RangeAddRangeMaxWithIndex(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_max_with_index<Monoid>, seg_act_add_with_index<Monoid, Action>, seg_comp_add<Action>,
        make_pair(-numeric_limits<Monoid>::max()/2, -1), -1);
};
template<class Monoid, class Action> LazySegmentTree<pair<Monoid,long long>, Action> RangeChangeRangeSum(int N = 0) {
    return LazySegmentTree<Monoid, Action>(
        N, seg_op_add_with_sum<Monoid>, seg_act_change_with_sum<Monoid, Action>, seg_comp_change<Action>,
        make_pair(Monoid(0), 1), -1);
};


//------------------------------//
// Tree
//------------------------------//

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


//------------------------------//
// Solver
//------------------------------//

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    // FastRead Read;
    // FastWrite Write;
    
    
}