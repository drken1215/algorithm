//
// Euclid 環上の行列 (加法・減法・乗法, 行列累乗, 行列式 (in O(N^3 log M))
//   Ring は「加法」「減法」「乗法」が定義されているクラスで、「除法」もできる（商が求められるもの）。コンストラクタで以下の情報を渡す。
//　　・ADD (加法), SUB (減法), MUL (乗法), DIV (商を求める), ADD_IDENTITY (加法の単位元), MUL_IDENTITY (乗法の単位元)
//   行列式を除算なしで O(N^4) で求める
//
// 例：
// 　 ・任意 mod の MintMatrx の行列式
// 　 ・mod 998244353 係数多項式の行列の行列式
//
// reference:
//   https://noshi91.hatenablog.com/entry/2020/11/28/115621
//
// verified:
//   Yosupo Library Checker - Determinant of Matrix (Arbitrary Mod)
//     https://judge.yosupo.jp/problem/matrix_det_arbitrary_mod
//
//   ABC 412 G - Degree Harmony（ただし、TLE）
//     https://atcoder.jp/contests/abc412/tasks/abc412_g
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Utility
//------------------------------//

// min non-negative i such that (x & (1 << i)) != 0
int bsf(int x) { return __builtin_ctz(x); }
int bsf(unsigned int x) { return __builtin_ctz(x); }
int bsf(long long x) { return __builtin_ctzll(x); }
int bsf(unsigned long long x) { return __builtin_ctzll(x); }

// min non-negative i such that n <= 2^i
template<class T> T ceil_pow2(T n) {
    T i = 0;
    while ((T(1) << i) < T(n)) i++;
    return i;
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


//------------------------------//
// Matrix
//------------------------------//

// general Euclid ring matrix (define ADD, SUB, MUL, DIV, ADD_IDENTITY, MUL_IDENTITY)
template<class Ring> struct EuclidRingMatrix {
    using FuncOperator = function<Ring(Ring, Ring)>;

    // inner value
    int H, W;
    vector<vector<Ring>> val;

    // operators
    FuncOperator ADD, SUB, MUL, DIV;
    Ring ADD_IDENTITY, MUL_IDENTITY;
    
    // constructors
    EuclidRingMatrix() {}
    EuclidRingMatrix(const EuclidRingMatrix&) = default;
    EuclidRingMatrix& operator = (const EuclidRingMatrix&) = default;
    EuclidRingMatrix(int h, int w
    , FuncOperator add, FuncOperator sub, FuncOperator mul, FuncOperator div
    , Ring add_id, Ring mul_id)
        : H(h), W(w), val(h, vector<Ring>(w, add_id))
        , ADD(add), SUB(sub), MUL(mul), DIV(div)
        , ADD_IDENTITY(add_id), MUL_IDENTITY(mul_id) {}
    void init(int h, int w
    , FuncOperator add, FuncOperator sub, FuncOperator mul, FuncOperator div
    , Ring add_id, Ring mul_id) {
        H = h, W = w;
        ADD = add, SUB = sub, MUL = mul, DIV = div;
        ADD_IDENTITY = add_id, MUL_IDENTITY = mul_id;
        val.assign(h, vector<Ring>(w, ADD_IDENTITY));
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
    vector<Ring>& operator [] (int i) { return val[i]; }
    const vector<Ring>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const EuclidRingMatrix &mat) {
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
    constexpr bool operator == (const EuclidRingMatrix &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const EuclidRingMatrix &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr EuclidRingMatrix& operator += (const EuclidRingMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = ADD(val[i][j], r.val[i][j]);
        return *this;
    }
    constexpr EuclidRingMatrix& operator -= (const EuclidRingMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = SUB(val[i][j], r.val[i][j]);
        return *this;
    }
    constexpr EuclidRingMatrix& operator *= (const Ring &v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = MUL(val[i][j], v);
        return *this;
    }
    constexpr EuclidRingMatrix& operator *= (const EuclidRingMatrix &r) {
        assert(width() == r.height());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        EuclidRingMatrix<Ring> res(height(), r.width(), ADD, SUB, MUL, DIV, ADD_IDENTITY, MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] = ADD(res[i][j], MUL(val[i][k], r.val[k][j]));
        return (*this) = res;
    }
    constexpr EuclidRingMatrix operator + () const { 
        return EuclidRingMatrix(*this);
    }
    constexpr EuclidRingMatrix operator + (const EuclidRingMatrix &r) const { 
        return EuclidRingMatrix(*this) += r;
    }
    constexpr EuclidRingMatrix operator - () const {
        EuclidRingMatrix res(*this);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res.val[i][j] = SUB(ADD_IDENTITY, res.val[i][j]);
        return res;
    }
    constexpr EuclidRingMatrix operator * (const Ring &v) const { 
        return EuclidRingMatrix(*this) *= v;
    }
    constexpr EuclidRingMatrix operator * (const EuclidRingMatrix &r) const { 
        return EuclidRingMatrix(*this) *= r;
    }
    constexpr vector<Ring> operator * (const vector<Ring> &v) const {
        assert(width() == v.size());
        vector<Ring> res(height(), ADD_IDENTITY);
        for (int i = 0; i < height(); i++)
            for (int j = 0; j < width(); j++)
                res[i] = ADD(res[i], MUL(val[i][j], v[j]));
        return res;
    }

    // transpose
    constexpr EuclidRingMatrix trans() const {
        EuclidRingMatrix<Ring> res(width(), height(), ADD, SUB, MUL, DIV, ADD_IDENTITY, MUL_IDENTITY);
        for (int row = 0; row < width(); row++)
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend constexpr EuclidRingMatrix trans(const EuclidRingMatrix &mat) {
        return mat.trans();
    }
    
    // pow
    constexpr EuclidRingMatrix pow(long long n) const {
        assert(height() == width());
        EuclidRingMatrix<Ring> res(height(), width(), ADD, SUB, MUL, DIV, ADD_IDENTITY, MUL_IDENTITY);
        EuclidRingMatrix<Ring> mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = MUL_IDENTITY;
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr EuclidRingMatrix pow(const EuclidRingMatrix &mat, long long n) {
        return mat.pow(n);
    }

    // determinant (without division, O(N^4))
    constexpr int find_pivot(int cur_rank, int col) const {
        int pivot = -1;
        for (int row = cur_rank; row < height(); ++row) {
            if (val[row][col] != ADD_IDENTITY) {
                pivot = row;
                break;
            }
        }
        return pivot;
    }
    constexpr Ring det() const {
        assert(height() == width());
        if (height() == 0) return MUL_IDENTITY;
        EuclidRingMatrix<Ring> A(*this);
        int rank = 0;
        Ring res = MUL_IDENTITY;
        for (int col = 0; col < width(); ++col) {
            int pivot = A.find_pivot(rank, col);
            if (pivot == -1) return ADD_IDENTITY;
            if (pivot != rank) swap(A[pivot], A[rank]), res = SUB(ADD_IDENTITY, res);
            for (int row = rank + 1; row < height(); ++row) {
                while (A[row][col] != ADD_IDENTITY) {
                    swap(A[rank], A[row]), res = SUB(ADD_IDENTITY, res);
                    Ring quo = DIV(A[row][col], A[rank][col]);
                    for (int col2 = rank; col2 < width(); ++col2) {
                        A[row][col2] = SUB(A[row][col2], MUL(A[rank][col2], quo));
                    }
                }
            }
            rank++;
        }
        for (int col = 0; col < height(); ++col) res = MUL(res, A[col][col]);
        return res;
    }
    friend constexpr Ring det(const EuclidRingMatrix &mat) {
        return mat.det();
    }
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


//------------------------------//
// NTT
//------------------------------//

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


//------------------------------//
// FPS
//------------------------------//

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


//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Determinant of Matrix (Arbitrary Mod)
void Yosupo_Determinant_of_Matrix_Arbitrary_Mod() {
    using mint = DynamicModint;
    int N, mod;
    cin >> N >> mod;
    mint::set_mod(mod);

    auto add = [&](mint a, mint b) -> mint { return a + b; };
    auto sub = [&](mint a, mint b) -> mint { return a - b; };
    auto mul = [&](mint a, mint b) -> mint { return a * b; };
    auto div = [&](mint a, mint b) -> mint { return (long long)(a.get()) / b.get(); };
    EuclidRingMatrix<mint> A(N, N, add, sub, mul, div, 0, 1);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) cin >> A[i][j];
    auto res = det(A);
    cout << res << endl;
}

// ABC 412 G - Degree Harmony
void ABC_412_G() {
    int N, M;
    cin >> N >> M;
    vector<int> A(N), id(N+1, 0), u(M), v(M);
    for (int i = 0; i < N; i++) cin >> A[i], id[i+1] = id[i] + A[i];
    for (int i = 0; i < M; i++) cin >> u[i] >> v[i], u[i]--, v[i]--;
    int V = id.back();

    const int MOD = 998244353;
    using mint = Fp<MOD>;
    auto add = [&](const FPS<mint> &a, const FPS<mint> &b) -> FPS<mint> { return a + b; };
    auto sub = [&](const FPS<mint> &a, const FPS<mint> &b) -> FPS<mint> { return a - b; };
    auto mul = [&](const FPS<mint> &a, const FPS<mint> &b) -> FPS<mint> { return a * b; };
    auto div = [&](const FPS<mint> &a, const FPS<mint> &b) -> FPS<mint> { return a / b; };
    EuclidRingMatrix<FPS<mint>> L(V, V, add, sub, mul, div, FPS<mint>(), FPS<mint>{1});
    for (int i = 0; i < N; i++) {
        for (int j = id[i]; j < id[i+1]; j++) {
            for (int k = id[i]; k < id[i+1]; k++) {
                if (j >= k) continue;
                int val = rand_int(0, MOD-1);
                L[j][k] = FPS<mint>{val};
                L[k][j] = FPS<mint>{-val};
            }
        }
    }
    for (int i = 0; i < M; i++) {
        for (int j = id[u[i]]; j < id[u[i]+1]; j++) {
            for (int k = id[v[i]]; k < id[v[i]+1]; k++) {
                int val = rand_int(0, MOD-1);
                L[j][k] = FPS<mint>{0, val};
                L[k][j] = FPS<mint>{0, -val};
            }
        }
    }
    auto f = det(L);
    int res = -1;
    for (int i = 0; i < f.size(); i++) if (f[i] != mint(0)) {
        res = i/2;
        break;
    }
    cout << res << endl;
}


int main() {
    //Yosupo_Determinant_of_Matrix_Arbitrary_Mod();
    ABC_412_G();
}