//
// 一般の可換環上の行列 (加法・乗法, 行列累乗, 行列式 (in O(N^4))
//   Ring は「加法」「減法」「乗法」が定義されているクラス。コンストラクタで以下の情報を渡す。
//　　・ADD (加法), SUB (減法), MUL (乗法), ADD_IDENTITY (加法の単位元), MUL_IDENTITY (乗法の単位元)
//　　・何もしなければ、通常の演算子「+」「-」「*」が呼び出される
//   行列式を除算なしで O(N^4) で求める
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


// general semiring matrix (define ADD, SUB, MUL, ADD_IDENTITY, MUL_IDENTITY)
template<class Ring> struct RingMatrix {
    using FuncOperator = function<Ring(Ring, Ring)>;

    // inner value
    int H, W;
    vector<vector<Ring>> val;

    // operators
    FuncOperator ADD, SUB, MUL;
    Ring ADD_IDENTITY;
    Ring MUL_IDENTITY;
    
    // constructors
    RingMatrix() {}
    RingMatrix(int h, int w
    , FuncOperator add, FuncOperator sub, FuncOperator mul, Ring add_id, Ring mul_id)
        : H(h), W(w), val(h, vector<Ring>(w, add_id))
        , ADD(add), SUB(sub), MUL(mul)
        , ADD_IDENTITY(add_id), MUL_IDENTITY(mul_id) {}
    void init(int h, int w
    , FuncOperator add, FuncOperator sub, FuncOperator mul, Ring add_id, Ring mul_id) {
        H = h, W = w;
        ADD = add, SUB = sub, MUL = mul;
        ADD_IDENTITY = add_id, MUL_IDENTITY = mul_id;
        val.assign(h, vector<Ring>(w, ADD_IDENTITY));
    }
    void resize(int h, int w) {
        H = h, W = w;
        val.resize(h);
        for (int i = 0; i < h; ++i) val[i].resize(w);
    }
    RingMatrix(const RingMatrix&) = default;
    RingMatrix& operator = (const RingMatrix&) = default;
    
    // getter and debugger
    constexpr int height() const { return H; }
    constexpr int width() const { return W; }
    constexpr bool empty() const { return height() == 0; }
    vector<Ring>& operator [] (int i) { return val[i]; }
    const vector<Ring>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const RingMatrix &mat) {
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
    constexpr bool operator == (const RingMatrix &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const RingMatrix &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr RingMatrix& operator += (const RingMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = ADD(val[i][j], r.val[i][j]);
        return *this;
    }
    constexpr RingMatrix& operator -= (const RingMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = SUB(val[i][j], r.val[i][j]);
        return *this;
    }
    constexpr RingMatrix& operator *= (const Ring &v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = MUL(val[i][j], v);
        return *this;
    }
    constexpr RingMatrix& operator *= (const RingMatrix &r) {
        assert(width() == r.height());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        RingMatrix<Ring> res(height(), r.width(), ADD, SUB, MUL, ADD_IDENTITY, MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] = ADD(res[i][j], MUL(val[i][k], r.val[k][j]));
        return (*this) = res;
    }
    constexpr RingMatrix operator + () const { 
        return RingMatrix(*this);
    }
    constexpr RingMatrix operator + (const RingMatrix &r) const { 
        return RingMatrix(*this) += r;
    }
    constexpr RingMatrix operator - () const {
        RingMatrix res(*this);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res.val[i][j] = SUB(ADD_IDENTITY, res.val[i][j]);
        return res;
    }
    constexpr RingMatrix operator * (const Ring &v) const { 
        return RingMatrix(*this) *= v;
    }
    constexpr RingMatrix operator * (const RingMatrix &r) const { 
        return RingMatrix(*this) *= r;
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
    constexpr RingMatrix trans() const {
        RingMatrix<Ring> res(width(), height(), ADD, SUB, MUL, ADD_IDENTITY, MUL_IDENTITY);
        for (int row = 0; row < width(); row++)
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend constexpr RingMatrix trans(const RingMatrix &mat) {
        return mat.trans();
    }
    
    // pow
    constexpr RingMatrix pow(long long n) const {
        assert(height() == width());
        RingMatrix<Ring> res(height(), width(), ADD, SUB, MUL, ADD_IDENTITY, MUL_IDENTITY);
        RingMatrix<Ring> mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = MUL_IDENTITY;
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr RingMatrix pow(const RingMatrix &mat, long long n) {
        return mat.pow(n);
    }

    // determinant (without division, O(N^4))
    constexpr Ring det() const {
        assert(height() == width());
        if (height() == 0) return MUL_IDENTITY;
        int N = height();
        vector<vector<Ring>> dp(N + 1, vector<Ring>(N + 1, ADD_IDENTITY));
        for (int i = 0; i <= N; i++) dp[i][i] = MUL_IDENTITY;
        for (int step = 0; step < N; step++) {
            vector<vector<Ring>> nex(N + 1, vector<Ring>(N + 1, ADD_IDENTITY));
            for (int row = 0; row < N; row++) {
                for (int col = row; col < N; col++) {
                    for (int col2 = row + 1; col2 < N; col2++) {
                        nex[row][col2] = SUB(nex[row][col2], MUL(dp[row][col], (*this)[col][col2]));
                    }
                    Ring tmp = MUL(dp[row][col], (*this)[col][row]);
                    for (int col2 = row + 1; col2 <= N; col2++) {
                        nex[col2][col2] = ADD(nex[col2][col2], tmp);
                    }
                }
            }
            swap(dp, nex);
        }
        return dp[N][N];
    }
    friend constexpr Ring det(const RingMatrix &mat) {
        return mat.det();
    }
};


//------------------------------//
// Examples
//------------------------------//

// yukicoder No.1303 Inconvenient Kingdom
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
};

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

void yukicoder_1303_general_det() {
    using mint = Fp<>;
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
        auto add = [&](mint a, mint b) -> mint { return a + b; };
        auto sub = [&](mint a, mint b) -> mint { return a - b; };
        auto mul = [&](mint a, mint b) -> mint { return a * b; };
        vector<int> conv(N, -1);
        int iter = 0;
        for (auto v : group) conv[v] = iter++;
        RingMatrix<mint> L(iter, iter, add, sub, mul, 0, 1);
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
        using Node = pair<mint, mint>;
        auto add = [&](Node a, Node b) -> Node { 
            return Node(a.first + b.first, a.second + b.second); 
        };
        auto sub = [&](Node a, Node b) -> Node { 
            return Node(a.first - b.first, a.second - b.second); 
        };
        auto mul = [&](Node a, Node b) -> Node { 
            return Node(a.first * b.first, a.first * b.second + a.second * b.first);
        };
        long long huben = 0;
        Node zero(0, 0), one(1, 0);
        RingMatrix<Node> L(N, N, add, sub, mul, zero, one);
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