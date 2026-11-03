//
// 一般の可換体上の行列 (加法・乗法, 行列累乗, 行列式 (in O(N^3))
//   Ring は「加法」「減法」「乗法」「除法」が定義されているクラス。コンストラクタで以下の情報を渡す。
//　　・ADD (加法), SUB (減法), MUL (乗法), DIV (除法), ADD_IDENTITY (加法の単位元), MUL_IDENTITY (乗法の単位元)
//
// verified:
//   yukicoder No.1303 Inconvenient Kingdom
//     https://yukicoder.me/problems/no/1303
//


#include <bits/stdc++.h>
using namespace std;


// general ring matrix (define ADD, SUB, MUL, DIV, ADD_IDENTITY, MUL_IDENTITY)
template<class Field> struct FieldMatrix {
    using FuncOperator = function<Field(Field, Field)>;

    // inner value
    int H, W;
    vector<vector<Field>> val;

    // operators
    FuncOperator ADD, SUB, MUL, DIV;
    Field ADD_IDENTITY, MUL_IDENTITY;
    
    // constructors
    FieldMatrix() {}
    FieldMatrix(const FieldMatrix&) = default;
    FieldMatrix& operator = (const FieldMatrix&) = default;
    FieldMatrix(int h, int w
    , FuncOperator add, FuncOperator sub, FuncOperator mul, FuncOperator div
    , Field add_id, Field mul_id)
        : H(h), W(w), val(h, vector<Field>(w, add_id))
        , ADD(add), SUB(sub), MUL(mul), DIV(div)
        , ADD_IDENTITY(add_id), MUL_IDENTITY(mul_id) {}
    void init(int h, int w
    , FuncOperator add, FuncOperator sub, FuncOperator mul, FuncOperator div
    , Field add_id, Field mul_id) {
        H = h, W = w;
        ADD = add, SUB = sub, MUL = mul, DIV = div;
        ADD_IDENTITY = add_id, MUL_IDENTITY = mul_id;
        val.assign(h, vector<Field>(w, ADD_IDENTITY));
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
    vector<Field>& operator [] (int i) { return val[i]; }
    const vector<Field>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const FieldMatrix &mat) {
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
    constexpr bool operator == (const FieldMatrix &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const FieldMatrix &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr FieldMatrix& operator += (const FieldMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = ADD(val[i][j], r.val[i][j]);
        return *this;
    }
    constexpr FieldMatrix& operator -= (const FieldMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = SUB(val[i][j], r.val[i][j]);
        return *this;
    }
    constexpr FieldMatrix& operator *= (const Field &v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = MUL(val[i][j], v);
        return *this;
    }
    constexpr FieldMatrix& operator *= (const FieldMatrix &r) {
        assert(width() == r.height());
        assert(ADD_IDENTITY == r.ADD_IDENTITY);
        assert(MUL_IDENTITY == r.MUL_IDENTITY);
        FieldMatrix<Field> res(height(), r.width(), ADD, SUB, MUL, DIV, ADD_IDENTITY, MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] = ADD(res[i][j], MUL(val[i][k], r.val[k][j]));
        return (*this) = res;
    }
    constexpr FieldMatrix operator + () const { 
        return FieldMatrix(*this);
    }
    constexpr FieldMatrix operator + (const FieldMatrix &r) const { 
        return FieldMatrix(*this) += r;
    }
    constexpr FieldMatrix operator - () const {
        FieldMatrix res(*this);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res.val[i][j] = SUB(ADD_IDENTITY, res.val[i][j]);
        return res;
    }
    constexpr FieldMatrix operator * (const Field &v) const { 
        return FieldMatrix(*this) *= v;
    }
    constexpr FieldMatrix operator * (const FieldMatrix &r) const { 
        return FieldMatrix(*this) *= r;
    }
    constexpr vector<Field> operator * (const vector<Field> &v) const {
        assert(width() == v.size());
        vector<Field> res(height(), ADD_IDENTITY);
        for (int i = 0; i < height(); i++)
            for (int j = 0; j < width(); j++)
                res[i] = ADD(res[i], MUL(val[i][j], v[j]));
        return res;
    }

    // transpose
    constexpr FieldMatrix trans() const {
        FieldMatrix<Field> res(width(), height(), ADD, SUB, MUL, ADD_IDENTITY, MUL_IDENTITY);
        for (int row = 0; row < width(); row++)
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend constexpr FieldMatrix trans(const FieldMatrix &mat) {
        return mat.trans();
    }
    
    // pow
    constexpr FieldMatrix pow(long long n) const {
        assert(height() == width());
        FieldMatrix<Field> res(height(), width(), ADD, SUB, MUL, ADD_IDENTITY, MUL_IDENTITY);
        FieldMatrix<Field> mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = MUL_IDENTITY;
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr FieldMatrix pow(const FieldMatrix &mat, long long n) {
        return mat.pow(n);
    }

    // gauss-jordan
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
    constexpr void sweep(int cur_rank, int col, int pivot, bool sweep_upper = true) {
        swap(val[pivot], val[cur_rank]);
        auto ifac = DIV(MUL_IDENTITY, val[cur_rank][col]);
        for (int col2 = cur_rank; col2 < width(); ++col2) {
            val[cur_rank][col2] = MUL(val[cur_rank][col2], ifac);
        }
        int row_start = (sweep_upper ? 0 : cur_rank + 1);
        for (int row = row_start; row < height(); ++row) {
            if (row != cur_rank && val[row][col] != ADD_IDENTITY) {
                auto fac = val[row][col];
                for (int col2 = cur_rank; col2 < width(); ++col2) {
                    val[row][col2] = SUB(val[row][col2], MUL(val[cur_rank][col2], fac));
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
    friend constexpr int gauss_jordan(FieldMatrix &mat, int not_sweep_width = 0, bool sweep_upper = true) {
        return mat.gauss_jordan(not_sweep_width, sweep_upper);
    }

    // rank
    constexpr int get_rank() const {
        if (height() == 0 || width() == 0) return 0;
        FieldMatrix A(*this);
        if (height() < width()) A = A.trans();
        return A.gauss_jordan(0, false);
    }
    friend constexpr int get_rank(const FieldMatrix &mat) {
        return mat.get_rank();
    }

    // find one solution
    friend constexpr int linear_equation
    (const FieldMatrix &mat, const vector<Field> &b, vector<Field> &res) {
        // extend
        FieldMatrix<Field> A(mat.height(), mat.width() + 1
        , mat.ADD, mat.SUB, mat.MUL, mat.DIV, mat.ADD_IDENTITY, mat.MUL_IDENTITY);
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
    friend constexpr int linear_equation(const FieldMatrix &mat, const vector<Field> &b) {
        vector<Field> res;
        return linear_equation(mat, b, res);
    }

    // find all solutions
    friend int linear_equation
    (const FieldMatrix &mat, const vector<Field> &b, vector<Field> &res, vector<vector<Field>> &zeros) {
        // extend
        FieldMatrix<Field> A(mat.height(), mat.width() + 1
        , mat.ADD, mat.SUB, mat.MUL, mat.DIV, mat.ADD_IDENTITY, mat.MUL_IDENTITY);
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) A[i][j] = mat.val[i][j];
            A[i].back() = b[i];
        }
        vector<int> core;
        int rank = A.gauss_jordan(core, 1);
        
        // check if it has no solution
        for (int row = rank; row < mat.height(); ++row) {
            if (A[row].back() != mat.ADD_IDENTITY) return -1;
        }

        // construct the core solution
        res.assign(mat.width(), mat.ADD_IDENTITY);
        for (int i = 0; i < (int)core.size(); i++) res[core[i]] = A[i].back();
    
        // construct the all solutions
        zeros.clear();
        vector<bool> use(mat.width(), 0);
        for (auto c : core) use[c] = true;
        for (int j = 0; j < mat.width(); j++) {
            if (use[j]) continue;
            vector<Field> zero(mat.width(), mat.ADD_IDENTITY);
            zero[j] = mat.MUL_IDENTITY;
            for (int i = 0; i < (int)core.size(); i++) zero[core[i]] = mat.SUB(mat.ADD_IDENTITY, A[i][j]);
            zeros.push_back(zero);
        }
        return rank;
    }
    
    // determinant
    constexpr Field det() const {
        assert(height() == width());
        if (height() == 0) return MUL_IDENTITY;
        FieldMatrix<Field> A(*this);
        int rank = 0;
        Field res = MUL_IDENTITY;
        for (int col = 0; col < width(); ++col) {
            int pivot = A.find_pivot(rank, col);
            if (pivot == -1) return ADD_IDENTITY;
            if (pivot != rank) res = SUB(ADD_IDENTITY, res);
            res = MUL(res, A[pivot][rank]);
            A.sweep(rank++, col, pivot, false);
        }
        return res;
    }
    friend constexpr Field det(const FieldMatrix &mat) {
        return mat.det();
    }

    // inv
    constexpr FieldMatrix inv() const {
        assert(height() == width());

        // extend
        FieldMatrix<Field> A(height(), width() + height());
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) A[i][j] = val[i][j];
            A[i][i+width()] = MUL_IDENTITY;
        }
        vector<int> core;
        int rank = A.gauss_jordan(height(), true);

        // gauss jordan
        if (rank < height()) return FieldMatrix();
        FieldMatrix<Field> res(height(), width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res[i][j] = A[i][j+width()];
        return res;
    }
    friend constexpr FieldMatrix inv(const FieldMatrix &mat) {
        return mat.inv();
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
        auto div = [&](mint a, mint b) -> mint { return a / b; };
        vector<int> conv(N, -1);
        int iter = 0;
        for (auto v : group) conv[v] = iter++;
        FieldMatrix<mint> L(iter, iter, add, sub, mul, div, 0, 1);
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
        auto div = [&](Node a, Node b) -> Node { 
            mint iv = b.first.inv();
            Node binv = Node(iv, -iv * iv * b.second);
            return mul(a, binv);
        };
        long long huben = 0;
        Node zero(0, 0), one(1, 0);
        FieldMatrix<Node> L(N, N, add, sub, mul, div, zero, one);
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