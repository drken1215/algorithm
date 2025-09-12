//
// 除算なし行列式, O(N^4)
//
// reference:
//   https://noshi91.hatenablog.com/entry/2020/11/28/115621
//
// verified:
//   yukicoder No.1303 Inconvenient Kingdom
//     https://yukicoder.me/problems/no/1303
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
// General Matrix
//------------------------------//

// modint matrix
template<class T> struct GeneralMatrix {
    // inner value
    int H, W;
    T ZERO = T();
    T ONE;
    vector<vector<T>> val;
    
    // constructors
    GeneralMatrix() : H(0), W(0) {}
    GeneralMatrix(int h, int w, const T &ZERO, const T &ONE)
        : H(h), W(w), ZERO(ZERO), ONE(ONE), val(h, vector<T>(w, ZERO)) {}
    GeneralMatrix(const GeneralMatrix &mat) 
        : H(mat.H), W(mat.W), ZERO(mat.ZERO), ONE(mat.ONE), val(mat.val) {}
    void init(int h, int w, const T &x) {
        H = h, W = w;
        val.assign(h, vector<T>(w, x));
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
    vector<T>& operator [] (int i) { return val[i]; }
    constexpr const vector<T>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const GeneralMatrix<T> &mat) {
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
    constexpr bool operator == (const GeneralMatrix &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const GeneralMatrix &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr GeneralMatrix& operator += (const GeneralMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] + r.val[i][j];
        return *this;
    }
    constexpr GeneralMatrix& operator -= (const GeneralMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] - r.val[i][j];
        return *this;
    }
    constexpr GeneralMatrix& operator *= (const T &v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] * v;
        return *this;
    }
    constexpr GeneralMatrix& operator *= (const GeneralMatrix &r) {
        assert(width() == r.height());
        GeneralMatrix<T> res(height(), r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] = res[i][j] + val[i][k] * r.val[k][j];
        return (*this) = res;
    }
    constexpr GeneralMatrix operator + () const { return GeneralMatrix(*this); }
    constexpr GeneralMatrix operator - () const { return GeneralMatrix(*this) = -GeneralMatrix(*this); }
    constexpr GeneralMatrix operator + (const GeneralMatrix &r) const { return GeneralMatrix(*this) += r; }
    constexpr GeneralMatrix operator - (const GeneralMatrix &r) const { return GeneralMatrix(*this) -= r; }
    constexpr GeneralMatrix operator * (const T &v) const { return GeneralMatrix(*this) *= v; }
    constexpr GeneralMatrix operator * (const GeneralMatrix &r) const { return GeneralMatrix(*this) *= r; }
    constexpr vector<T> operator * (const vector<T> &v) const {
        assert(width() == v.size());
        vector<T> res(height(), ZERO);
        for (int i = 0; i < height(); i++)
            for (int j = 0; j < width(); j++)
                res[i] += val[i][j] * v[j];
        return res;
    }

    // transpose
    constexpr GeneralMatrix trans() const {
        GeneralMatrix<T> res(width(), height());
        for (int row = 0; row < width(); row++) 
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend constexpr GeneralMatrix<T> trans(const GeneralMatrix<T> &mat) {
        return mat.trans();
    }
    
    // pow
    constexpr GeneralMatrix pow(long long n) const {
        assert(height() == width());
        GeneralMatrix<T> res(height(), width()), mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = ONE;
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr GeneralMatrix<T> pow(const GeneralMatrix<T> &mat, long long n) {
        return mat.pow(n);
    }
    
    // determinant (without division, O(N^4))
    constexpr T det() const {
        assert(height() == width());
        if (height() == 0) return ONE;
        int N = height();
        vector<vector<T>> dp(N + 1, vector<T>(N + 1, ZERO));
        for (int i = 0; i <= N; i++) dp[i][i] = ONE;
        for (int step = 0; step < N; step++) {
            vector<vector<T>> nex(N + 1, vector<T>(N + 1, ZERO));
            for (int row = 0; row < N; row++) {
                for (int col = row; col < N; col++) {
                    for (int col2 = row + 1; col2 < N; col2++) {
                        nex[row][col2] = nex[row][col2] - dp[row][col] * (*this)[col][col2];
                    }
                    T tmp = dp[row][col] * (*this)[col][row];
                    for (int col2 = row + 1; col2 <= N; col2++) {
                        nex[col2][col2] = nex[col2][col2] + tmp;
                    }
                }
            }
            swap(dp, nex);
        }
        return dp[N][N];
    }
    friend constexpr T det(const GeneralMatrix<T> &mat) {
        return mat.det();
    }

    // determinant (by Euclidean Algorithm)
    constexpr int find_pivot(int cur_rank, int col) const {
        int pivot = -1;
        for (int row = cur_rank; row < height(); ++row) {
            if (val[row][col] != ZERO) {
                pivot = row;
                break;
            }
        }
        return pivot;
    }
    constexpr T det_euclid() const {
        assert(height() == width());
        if (height() == 0) return ONE;
        GeneralMatrix<T> A(*this);
        int rank = 0;
        T res = ONE;
        for (int col = 0; col < width(); ++col) {
            int pivot = A.find_pivot(rank, col);
            if (pivot == -1) return ZERO;
            if (pivot != rank) swap(A[pivot], A[rank]), res = -res;
            for (int row = rank + 1; row < height(); ++row) {
                while (A[row][col] != ZERO) {
                    swap(A[rank], A[row]), res = -res;
                    auto quo = A[row][col] / A[rank][col];
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
    friend constexpr T det_euclid(const GeneralMatrix<T> &mat) {
        return mat.det_euclid();
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


//------------------------------//
// Examples
//------------------------------//

// yukicoder No.1303 Inconvenient Kingdom
using mint = Fp<>;
using Node = pair<mint, mint>;
Node operator + (Node a, Node b) { return Node(a.first + b.first, a.second + b.second); }
Node operator - (Node a, Node b) { return Node(a.first - b.first, a.second - b.second); }
Node operator * (Node a, Node b) { return Node(a.first * b.first, a.first * b.second + a.second * b.first); }
void yukicoder_1303_general_det() {
    int N, M, u, v;
    cin >> N >> M;
    vector G(N, vector(N, 0));
    vector degs(N, 0);
    UnionFind uf(N);
    for (int i = 0; i < M; i++) {
        cin >> u >> v, u--, v--;
        G[u][v]++, G[v][u]++, degs[u]++, degs[v]++;
        uf.merge(u, v);
    }

    auto calc = [&](const vector<int> &group) -> mint {
        vector<int> conv(N, -1);
        int iter = 0;
        for (auto v : group) conv[v] = iter++;
        GeneralMatrix<mint> L(iter, iter, mint(0), mint(1));
        for (int i = 0; i < iter; i++) L[i][i] = 1;
        for (auto v1 : group) {
            int i = conv[v1];
            int deg = 0;
            for (auto v2 : group) {
                if (G[v1][v2]) deg += G[v1][v2];
                int j = conv[v2];
                if (i < iter - 1 && j < iter - 1) L[i][j] = -G[v1][v2];
            }
            if (i < iter - 1) L[i][i] = deg;
        }
        return det(L);
    };

    auto groups = uf.groups();
    if (groups.size() > 1) {
        mint res = 1;
        vector<long long> siz;
        for (auto group : groups) siz.push_back(group.size()), res *= calc(group);
        sort(siz.begin(), siz.end(), greater<long long>());
        long long sum = siz[0] + siz[1], sum2 = sum * sum;
        for (int i = 2; i < (int)siz.size(); i++) sum += siz[i], sum2 += siz[i] * siz[i];
        long long huben = (sum * sum - sum2);
        if (siz[0] == siz[1]) {
            long long sumsiz = 0, sumsiz2 = 0;
            for (auto s : siz) if (s == siz[0]) sumsiz += s, sumsiz2 += s * s;
            long long fac = (sumsiz * sumsiz - sumsiz2) / 2;
            res *= fac;
        } else {
            long long sum_sub = 0;
            for (auto s : siz) if (s == siz[1]) sum_sub += s;
            long long fac = siz[0] * sum_sub;
            res *= fac;
        }
        cout << huben << endl << res << endl;
    } else {
        long long huben = 0;
        Node zero(0, 0), one(1, 0);
        GeneralMatrix<Node> L(N, N, zero, one);
        for (int i = 0; i < N; i++) L[i][i] = one;
        for (int i = 0; i < N-1; i++) {
            L[i][i] = Node(mint(degs[i]), mint(N - 1 - degs[i]));
            for (int j = 0; j < N-1; j++) {
                if (i == j) continue;
                if (G[i][j]) L[i][j] = Node(mint(-G[i][j]), 0);
                else L[i][j] = Node(0, mint(-1));
            }
        }
        Node f = det(L);
        mint res = f.first + f.second;
        cout << huben << endl << res << endl;
    }
}


int main() {
    yukicoder_1303_general_det();
}