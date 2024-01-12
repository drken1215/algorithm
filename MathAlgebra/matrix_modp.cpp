//
// mod. p 行列 (行列累乗、掃き出し法)
//
// verified:
//   TCO 2013 Round 2A Med TheMagicMatrix
//     https://vjudge.net/problem/TopCoder-12495
//
//   AOJ 3369 (?) Namori Counting (OUPC 2023 day2-D)
//     https://onlinejudge.u-aizu.ac.jp/beta/room.html#OUPC2023Day2/problems/D
//
//   AOJ 3369 (?) Namori Counting (OUPC 2023 day2-D)
//     https://onlinejudge.u-aizu.ac.jp/beta/room.html#OUPC2023Day2/problems/D
//


#include <bits/stdc++.h>
using namespace std;


// matrix
template<class mint> struct MintMatrix {
    // inner value
    vector<vector<mint>> val;
    
    // constructors
    MintMatrix(int H, int W, mint x = 0) : val(H, vector<mint>(W, x)) {}
    MintMatrix(const MintMatrix &mat) : val(mat.val) {}
    void init(int H, int W, mint x = 0) {
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
                val[i][j] += r[i][j];
            }
        }
        return *this;
    }
    constexpr MintMatrix& operator -= (const MintMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                val[i][j] -= r[i][j];
            }
        }
        return *this;
    }
    constexpr MintMatrix& operator *= (const mint &v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] *= v;
        return *this;
    }
    constexpr MintMatrix& operator *= (const MintMatrix &r) {
        assert(width() == r.height());
        MintMatrix<mint> res(height(), r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] += val[i][k] * r[k][j];
        return (*this) = res;
    }
    constexpr MintMatrix operator + () const { return MintMatrix(*this); }
    constexpr MintMatrix operator - () const { return MintMatrix(*this) *= mint(-1); }
    constexpr MintMatrix operator + (const MintMatrix &r) const { return MintMatrix(*this) += r; }
    constexpr MintMatrix operator - (const MintMatrix &r) const { return MintMatrix(*this) -= r; }
    constexpr MintMatrix operator * (const mint &v) const { return MintMatrix(*this) *= v; }
    constexpr MintMatrix operator * (const MintMatrix &r) const { return MintMatrix(*this) *= r; }
    
    // pow
    constexpr MintMatrix pow(long long n) const {
        assert(height() == width());
        MintMatrix<mint> res(height(), width()),  mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
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
    friend constexpr int gauss_jordan(MintMatrix<mint> &mat, int not_sweep_width = 0) {
        return mat.gauss_jordan(not_sweep_width);
    }
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

// modint
template<int MOD> struct Fp {
    // inner value
    long long val;
    
    // constructor
    constexpr Fp() : val(0) { }
    constexpr Fp(long long v) : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr Fp(const Fp &v) : val(v.get()) { }
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



//------------------------------//
// Examples
//------------------------------//

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


// TCO 2013 Round 2A Med TheMagicMatrix
class TheMagicMatrix {
public:
    int find(int n, vector<int> rows, vector<int> cols, vector<int> vals) {
        const int MOD = 1234567891;
        using mint = Fp<MOD>;
        using mint2 = Fp<2>;
        using mint5 = Fp<5>;
        
        // 数字のない行と列がある場合
        int m = (int)rows.size();
        set<int> sr, sc;
        for (int i = 0; i < m; ++i) {
            sr.insert(rows[i]);
            sc.insert(cols[i]);
        }
        if (sr.size() < n && sc.size() < n) {
            long long ex = (n - 1) * (n - 1) - m + 1;
            return mint(10).pow(ex).get();
        }

        // 連立一次方程式を立てる
        MintMatrix<mint2> M2(n * 2 + m, n * n);
        MintMatrix<mint5> M5(n * 2 + m, n * n);
        vector<mint2> b2(n * 2 + m);
        vector<mint5> b5(n * 2 + m);

        // 行和
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                int id = i * n + j;
                M2[i][id] = 1;
                M5[i][id] = 1;
            }
        }

        // 列和
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                int id = i * n + j;
                M2[j + n][id] = 1;
                M5[j + n][id] = 1;
            }
        }

        // 条件
        for (int k = 0; k < m; ++k) {
            int id = rows[k] * n + cols[k];
            M2[k + n * 2][id] = 1;
            M5[k + n * 2][id] = 1;
            b2[k + n * 2] = vals[k];
            b5[k + n * 2] = vals[k];
        }

        // X = 0, 1, ..., 9
        mint res = 0;
        for (int X = 0; X < 10; ++X) {
            for (int i = 0; i < n * 2; ++i) {
                b2[i] = X;
                b5[i] = X;
            }
            int rank2 = linear_equation(M2, b2);
            int rank5 = linear_equation(M5, b5);
            if (rank2 == -1 || rank5 == -1) continue;
            mint tmp = mint(2).pow(n * n - rank2) * mint(5).pow(n * n - rank5);
            res += tmp;
        }
        return res.get();
    }
};

void TCO_2013_Round2_A() {
    int n, m;
    cin >> n >> m;
    vector<int> r(m), c(m), v(m);
    for (int i = 0; i < m; ++i) cin >> r[i] >> c[i] >> v[i];
    TheMagicMatrix tmm;
    cout << tmm.find(n, r, c, v) << endl;
}


int main() {
    AOJ_3369();
    //TCO_2013_Round2_A();
}

