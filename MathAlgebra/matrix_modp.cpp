//
// mod. p 行列 (行列累乗、掃き出し法、逆行列、連立一次方程式)
//
// verified:
//   Yosupo Library Checker - Pow of Matrix
//     https://judge.yosupo.jp/problem/pow_of_matrix 
//
//   Yosupo Library Checker - Determinant of Matrix
//     https://judge.yosupo.jp/problem/matrix_det
//
//   Yosupo Library Checker - Determinant of Matrix (Arbitrary Mod)
//     https://judge.yosupo.jp/problem/matrix_det_arbitrary_mod
//
//   Yosupo Library Checker - Rank of Matrix
//     https://judge.yosupo.jp/problem/matrix_rank
//
//   Yosupo Library Checker - System of Linear Equations
//     https://judge.yosupo.jp/problem/system_of_linear_equations
//
//   Yosupo Library Checker - Inverse Matrix
//     https://judge.yosupo.jp/problem/inverse_matrix
//
//   AOJ 3369 (?) Namori Counting (OUPC 2023 day2-D)
//     https://onlinejudge.u-aizu.ac.jp/beta/room.html#OUPC2023Day2/problems/D
//
//   ARC 199 B - Adjacent Replace
//     https://atcoder.jp/contests/arc199/tasks/arc199_b
//


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


//------------------------------//
// Modint Matrix
//------------------------------//

// matrix
template<class mint> struct MintMatrix {
    // inner value
    int H, W;
    vector<vector<mint>> val;
    
    // constructors
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
};


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
// Examples
//------------------------------//

// Yosupo Library Checker - Pow of Matrix
void Yosupo_Pow_of_Matrix() {
    using mint = Fp<>;
    FastRead Read;
    FastWrite Write;
    long long N, K;
    Read(N, K);
    MintMatrix<mint> A(N, N);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Read(A[i][j].val);
    auto res = pow(A, K);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) Write(res[i][j].val), Write(' ');
        Write('\n');
    }
}

// Yosupo Library Checker - Determinant of Matrix
void Yosupo_Determinant_of_Matrix() {
    using mint = Fp<>;
    FastRead Read;
    FastWrite Write;
    int N;
    Read(N);
    MintMatrix<mint> A(N, N);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Read(A[i][j].val);
    auto res = det(A);
    Write(res.val), Write('\n');
}

// Yosupo Library Checker - Determinant of Matrix (Arbitrary Mod)
void Yosupo_Determinant_of_Matrix_Arbitrary_Mod() {
    using mint = DynamicModint;
    FastRead Read;
    FastWrite Write;
    int N, mod;
    Read(N, mod);
    mint::set_mod(mod);
    MintMatrix<mint> A(N, N);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Read(A[i][j].val);
    auto res = det_nonprime_mod(A);
    Write(res.val), Write('\n');
}

// Yosupo Library Checker - Rank of Matrix 
void Yosupo_Rank_of_Matrix() {
    using mint = Fp<>;
    FastRead Read;
    FastWrite Write;
    int N, M;
    Read(N, M);
    MintMatrix<mint> A(N, M);
    for (int i = 0; i < N; i++) for (int j = 0; j < M; j++) Read(A[i][j].val);
    auto res = get_rank(A);
    Write(res), Write('\n');
}

// Yosupo Library Checker - System of Linear Equations
void Yosupo_System_of_Linear_Equations() {

}

// Yosupo Library Checker - Inverse Matrix
void Yosupo_Inverse_Matrix() {

}

// AOJ 3369 Namori Counting (OUPC 2023 day2-D)
void AOJ_3369() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    int N, M;
    cin >> N >> M;
    vector<int> deg(N, 0);
    vector<vector<int>> G(N, vector<int>(N, 0));
    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        --u, --v;
        ++G[u][v], ++G[v][u];
        ++deg[u], ++deg[v];
    }
    
    // ラプラシアン行列の余因子を求めるため、行・列の末尾を削る
    MintMatrix<mint> L(N - 1, N - 1, 0);
    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            if (i == j) L[i][j] = deg[i];
            else L[i][j] = -G[i][j];
        }
    }
    mint res = det(L) * (M - N + 1);
    cout << res << endl;
}

//ARC 199 B - Adjacent Replace
void solve() {
    using mint = Fp<2>;
    auto construct = [&](vector<int> need) -> vector<int> {
        int N = need.size();
        vector<int> res;
        while (true) {
            bool finish = true;
            if (need[0] == 0) finish = false;
            for (int i = 1; i < N; i++) if (need[i] == 1) finish = false;
            if (finish) break;

            // (1, 1, ...) と (..., 1, 1) の解消
            if (need[0] == 1 && need[1] == 1) {
                if (need[2] <= 0) {
                    res.push_back(1);
                    need[0] = 1, need[1] = 0;
                } else {
                    res.push_back(2), res.push_back(1);
                    need[0] = 1, need[1] = 0, need[2] = 0;
                }
            }
            if (need[N-1] == 1 && need[N-2] == 1) {
                if (need[N-3] <= 0) {
                    res.push_back(N-1);
                    need[N-1] = 1, need[N-2] = 0;
                } else {
                    res.push_back(N-2), res.push_back(N-1);
                    need[N-1] = 1, need[N-2] = 0, need[N-3] = 0;
                }
            }

            // routine
            for (int i = 0; i+1 < N; i++) {
                if (need[i] == -1 && need[i+1] == -1) continue;
                if (need[i] == 1 && need[i+1] == 1) {
                    int j = i;
                    while (j < N && need[j] == need[i]) j++;
                    for (int k = j-2; k >= i; k--) {
                        res.push_back(k+1);
                        need[k+1] = 0;
                    }
                } else if (need[i] <= 0 && need[i+1] <= 0) {
                    res.push_back(i+1), res.push_back(i+1);
                    need[i] = need[i+1] = -1;
                } else if (need[i] == -1 && need[i+1] == 1) {
                    res.push_back(i+1);
                    need[i] = 1, need[i+1] = 0;
                } else if (i+2 < N && need[i] == 0 && need[i+1] == 1 && need[i+2] == -1) {
                    res.push_back(i+2), res.push_back(i+1), res.push_back(i+1), 
                    res.push_back(i+2), res.push_back(i+1);
                    need[i] = 1, need[i+1] = 0, need[i+2] = 0;
                }
            }
        }
        return res;
    };

    long long N, K;
    cin >> N >> K;
    vector<long long> A(N);
    MintMatrix<mint> M(60, N);
    vector<mint> v(60, 0);
    for (int i = 0; i < N; i++) {
        cin >> A[i];
        for (int d = 0; d < 60; d++) if (A[i] >> d & 1) M[d][i] = 1;
    }
    for (int d = 0; d < 60; d++) if (K >> d & 1) v[d] = 1;

    // find solutions
    vector<mint> ans;
    vector<vector<mint>> zeros;
    int rank = linear_equation(M, v, ans, zeros);
    if (rank == -1) {
        cout << "No" << endl;
        return;
    }
    bool exist = false;
    vector<int> need(N, 0);
    for (long long bit = 0; bit < (1LL<<zeros.size()); bit++) {
        vector<mint> tmp = ans;
        for (int i = 0; i < zeros.size(); i++) {
            if (bit >> i & 1) for (int j = 0; j < tmp.size(); j++) tmp[j] += zeros[i][j];
        }
        bool ok = false;
        for (int i = 0; i + 1 < tmp.size(); i++) if (tmp[i] == tmp[i+1]) ok = true;
        if (ok) {
            exist = true;
            for (int i = 0; i < N; i++) need[i] = tmp[i].val;
            break;
        }
    }
    if (!exist) {
        cout << "No" << endl;
        return;
    }

    auto check = [&](const vector<int> &res) -> bool {
        if (res.size() > 10000) return false;
        for (int i = 0; i < res.size(); i++) {
            int id = res[i];
            if (id < 1 || id > N - 1) return false;
            id--;
            long long x = A[id] ^ A[id + 1];
            A[id] = x, A[id + 1] = x;
        }
        return (A[0] == K);
    };

    auto res = construct(need);
    cout << "Yes" << endl;
    cout << res.size() << endl;
    for (auto val : res) cout << val << " ";
    cout << endl;
}
void ARC_199_B() {
    int T;
    cin >> T;
    while (T--) solve();
}


int main() {
    //Yosupo_Pow_of_Matrix();
    //Yosupo_Determinant_of_Matrix();
    //Yosupo_Determinant_of_Matrix_Arbitrary_Mod();
    Yosupo_Rank_of_Matrix();
    //Yosupo_System_of_Linear_Equations();
    //Yosupo_Inverse_Matrix();
    //AOJ_3369();
    //ARC_199_B();
}