//
// 一般の可換環上の行列 (加法・乗法, 行列累乗, 行列式 (in O(N^4))
//   Ring は「加法」「減法」「乗法」が定義されているクラス。コンストラクタで以下の情報を渡す。
//　　・コンストラクタで、ADD (加法), MUL (乗法), ADD_IDENTITY (加法の単位元), MUL_IDENTITY (乗法の単位元)
//　　・何もしなければ、通常の演算子「+」「-」「*」が呼び出される
//   行列式を除算なしで O(N^4) で求める方法がある。
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
    Ring ADD_IDENTITY = Ring(), MUL_IDENTITY = Ring(1);
    FuncOperator ADD = [](const Ring &a, const Ring &b) -> Ring { return a + b; };
    FuncOperator SUB = [](const Ring &a, const Ring &b) -> Ring { return a - b; };
    FuncOperator MUL = [](const Ring &a, const Ring &b) -> Ring { return a * b; };
    
    // constructors
    RingMatrix() : H(0), W(0) {}
    RingMatrix(int H, int W) : H(H), W(W), val(H, vector<Ring>(W, ADD_IDENTITY)) {}
    RingMatrix(int H, int W, Ring v) : H(H), W(W), val(H, vector<Ring>(W, v)) {}
    RingMatrix(const RingMatrix &mat) 
        : H(mat.H), W(mat.W), val(mat.val)
        , ADD(mat.ADD), SUB(mat.sub), MUL(mat.MUL)
        , ADD_IDENTITY(mat.ADD_IDENTITY), MUL_IDENTITY(mat.MUL_IDENTITY) {}
    RingMatrix(int H, int W
        , const FuncOperator add, const FuncOperator sub, const FuncOperator mul
        , const Ring &add_identity, const Ring &mul_identity) {
        init(H, W, add, mul, add_identity, mul_identity);
    }
    void init(int h, int w, const Ring &x) {
        H = h, W = w;
        val.assign(h, vector<Ring>(w, x));
    }
    void init(int h, int w
    , const FuncOperator add, const FuncOperator sub, const FuncOperator mul
    , const Ring &add_identity, const Ring &mul_identity) {
        H = h, W = w;
        ADD = add, SUB = sub, MUL = mul;
        ADD_IDENTITY = add_identity, MUL_IDENTITY = mul_identity;
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
    friend constexpr ostream& operator << (ostream &os, const RingMatrix<Ring> &mat) {
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
        assert(ADD_IDENTITY == r.ADD_IDENTITY), assert(MUL_IDENTITY == r.MUL_IDENTITY);
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                val[i][j] = ADD(val[i][j], r.val[i][j]);
            }
        }
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
        assert(ADD_IDENTITY == r.ADD_IDENTITY), assert(MUL_IDENTITY == r.MUL_IDENTITY);
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
    friend constexpr RingMatrix<Ring> trans(const RingMatrix<Ring> &mat) {
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
    friend constexpr RingMatrix<Ring> pow(const RingMatrix<Ring> &mat, long long n) {
        return mat.pow(n);
    }
};


//------------------------------//
// Examples
//------------------------------//

// yukicoder No.1303 Inconvenient Kingdom
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