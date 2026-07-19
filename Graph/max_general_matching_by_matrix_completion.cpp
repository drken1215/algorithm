//
// 一般グラフの最大マッチング (by 行列補完, in O(V^3))
//
// reference:
//   adamant: Randomized general matching with Tutte matrix
//     https://codeforces.com/blog/entry/92400
//
// verified:
//   Yosupo Library Checker - Matching on General Graph
//     https://judge.yosupo.jp/problem/general_matching
//
//   JOI春合宿2016 マッチングコンテスト
//     https://atcoder.jp/contests/joisc2016matching
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;

/*
    グラフ G に対する Tutte 行列 T
    ・辺 (i, j) があるとき
        T[i,j] = x[i,j] (i < j), -x[i,j] (i > j)
    ・辺 (i, j) がない or i = j のとき
        T[i,j] = 0

    このとき
    ・rank(T) = G の最大マッチングに含まれる頂点の個数（とくに、det(T) ≠ 0 　⇔ 　G が完全マッチングをもつ）
    ・T の行の部分集合がフルランク  ⇔  対応する G の部分グラフが完全マッチングをもつ
    ・G が完全マッチングをもつならば、「 T^{-1}[i,j] ≠ 0　　⇔ 　辺 (i, j) を含む完全マッチングが存在 」

    最大マッチングの復元方法
    1. Tutte 行列 T の独立な行集合を求め (以降、T をその行集合・列集合に制限した行列とする)、B = T^{-1} とする
    2. T[i,j] ≠ 0, B[i,j] ≠ 0 となる (i, j) を求め、それをマッチングに加える
    3. T, B から i, j 要素を削除する
    4. マッチングが求まるまで 2, 3 を繰り返す
*/


//------------------------------//
// Subroutine
//------------------------------//

// static modint
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

// modint matrix
template<class mint> struct MintMatrix {
    // inner value
    int H, W;
    vector<vector<mint>> val;
    
    // constructors
    MintMatrix() {}
    MintMatrix(const MintMatrix&) = default;
    MintMatrix& operator = (const MintMatrix&) = default;
    MintMatrix(int h, int w) : H(h), W(w), val(h, vector<mint>(w)) {}
    MintMatrix(int h, int w, mint x) : H(h), W(w), val(h, vector<mint>(w, x)) {}
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
    const vector<mint>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const MintMatrix &mat) {
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
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] + r.val[i][j];
        return *this;
    }
    constexpr MintMatrix& operator -= (const MintMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] - r.val[i][j];
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
    constexpr MintMatrix operator + () const { 
        return MintMatrix(*this);
    }
    constexpr MintMatrix operator - () const {
        MintMatrix res(*this);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res.val[i][j] = -res.val[i][j];
        return res;
    }
    constexpr MintMatrix operator + (const MintMatrix &r) const { 
        return MintMatrix(*this) += r;
    }
    constexpr MintMatrix operator - (const MintMatrix &r) const {
        return MintMatrix(*this) -= r;
    }
    constexpr MintMatrix operator * (const mint &v) const {
        return MintMatrix(*this) *= v;
    }
    constexpr MintMatrix operator * (const MintMatrix &r) const {
        return MintMatrix(*this) *= r;
    }
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
        for (int row = 0; row < width(); row++)
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend constexpr MintMatrix trans(const MintMatrix &mat) {
        return mat.trans();
    }
    
    // pow
    constexpr MintMatrix pow(long long n) const {
        assert(height() == width());
        MintMatrix<mint> res(height(), width());
        MintMatrix<mint> mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = mint(1);
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr MintMatrix pow(const MintMatrix &mat, long long n) {
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
    friend constexpr int gauss_jordan(MintMatrix &mat, int not_sweep_width = 0, bool sweep_upper = true) {
        return mat.gauss_jordan(not_sweep_width, sweep_upper);
    }

    // rank
    constexpr int get_rank() const {
        if (height() == 0 || width() == 0) return 0;
        MintMatrix A(*this);
        if (height() < width()) A = A.trans();
        return A.gauss_jordan(0, false);
    }
    friend constexpr int get_rank(const MintMatrix &mat) {
        return mat.get_rank();
    }

    // find one solution
    friend constexpr int linear_equation
    (const MintMatrix &mat, const vector<mint> &b, vector<mint> &res) {
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
    friend constexpr int linear_equation(const MintMatrix &mat, const vector<mint> &b) {
        vector<mint> res;
        return linear_equation(mat, b, res);
    }

    // find all solutions
    friend int linear_equation
    (const MintMatrix &mat, const vector<mint> &b, vector<mint> &res, vector<vector<mint>> &zeros) {
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
    friend constexpr mint det(const MintMatrix &mat) {
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
    friend constexpr mint det_nonprime_mod(const MintMatrix &mat) {
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
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res[i][j] = A[i][j+width()];
        return res;
    }
    friend constexpr MintMatrix inv(const MintMatrix &mat) {
        return mat.inv();
    }
};


//------------------------------//
// General Matching 
//------------------------------//

// find general matching of graph G
template<class mint> struct GeneralMatching {
    // rng
    constexpr static auto gm_rand_unsigned_int = []() -> unsigned int {
        static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
        unsigned int tt = (tx^(tx<<11));
        tx = ty; ty = tz; tz = tw;
        return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
    };
    constexpr static auto gm_rand_int = [](int minv, int maxv) -> int {
        return gm_rand_unsigned_int() % (maxv - minv + 1) + minv;
    };

    // inner values
    int MOD = mint::get_mod();
    int N;
    vector<vector<bool>> G;  // graph (adjacency list)

    // results
    MintMatrix<mint> T, iT;  // Tutte matrix, inv of Tutte matrix
    vector<int> bases;  // independent set

    // constructors
    GeneralMatching() {}
    GeneralMatching(const vector<vector<bool>> &graph) {
        init(graph);
    }
    constexpr void init(const vector<vector<bool>> &graph) {
        N = (int)graph.size();
        G = graph;
        T.init(N, N, mint(0)), iT.init(N, N, mint(0));
        for (int i = 0; i < N; i++) {
            assert(G[i].size() == N);
            iT[i][i] = mint(1);
            for (int j = 0; j < i; j++) {
                if (G[i][j]) {
                    auto x = gm_rand_int(1, MOD-1);
                    T[i][j] = x, T[j][i] = -x;
                }
            }
        }
    }

    // has perfect matching?
    constexpr bool has_perfect_matching() const {
        return det(T) != 0;
    }

    // find independent set of rows of T
    constexpr vector<int> find_independent_set() {
        bases.clear();
        vector<bool> used(N, false);
        for (int col = 0; col < N; col++) {
            int p = 0;
            while (p < N && (used[p] || T[p][col] == mint(0))) p++;
            if (p == N) continue;
            auto ifac = T[p][col].inv();
            bases.emplace_back(p);
            used[p] = true;
            for (int col2 = 0; col2 < N; col2++) {
                T[p][col2] *= ifac;
                iT[p][col2] *= ifac;
            }
            for (int row = 0; row < N; row++) {
                if (row != p && T[row][col] != mint(0)) {
                    auto fac = T[row][col];
                    for (int col2 = 0; col2 < N; col2++) {
                        T[row][col2] -= T[p][col2] * fac;
                        iT[row][col2] -= iT[p][col2] * fac;
                    }
                }
            }
        }
        sort(bases.begin(), bases.end());
        return bases;
    }

    // find max matching
    constexpr vector<pair<int,int>> solve() {
        auto T0 = T;
        const auto &bases = find_independent_set();

        // construct inverse submatrix
        int M = (int)bases.size();
        MintMatrix<mint> inv(M, M, mint(0));
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) {
            if (T[bases[j]][bases[i]] != mint(0)) {
                for (int k = 0; k < M; k++) {
                    inv[i][k] = iT[bases[j]][bases[k]];
                }
            }
        }

        // reconstruct max matching
        vector<pair<int,int>> res;
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) {
            int u = bases[i], v = bases[j];
            if (G[u][v] && inv[i][j] != mint(0)) {
                res.emplace_back(u, v);
                auto sij = inv[i][j].inv(), sji = inv[j][i].inv();
                for (int col = 0; col < M; col++) inv[i][col] *= sij; 
                for (int col = 0; col < M; col++) inv[j][col] *= sji;
                for (int k = 0; k < M; k++) {
                    if (k != i && k != j) {
                        auto fki = inv[k][i], fkj = inv[k][j];
                        for (int col = 0; col < M; col++) {
                            inv[k][col] -= fkj * inv[i][col] + fki * inv[j][col];
                        }
                    }
                }
                for (int col = 0; col < M; col++) inv[i][col] = inv[j][col] = mint(0);
            }
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Matching on General Graph
void yosupo_general_matching() {
    using mint = Fp<998244353>;
    int N, M, u, v;
    cin >> N >> M;
    vector<vector<bool>> G(N, vector<bool>(N, false));
    for (int i = 0; i < M; i++) {
        cin >> u >> v;
        G[u][v] = G[v][u] = true;
    }
    GeneralMatching<mint> gm(G);
    auto res = gm.solve();
    cout << res.size() << endl;
    for (auto [u, v] : res) cout << u << " " << v << endl;
}

// JOI春合宿2016 マッチングコンテスト
void JOI_2026_Matching() {
    using mint = Fp<>;
    int T;
    cin >> T;
    while (T--) {
        int N, M, a, b;
        cin >> N >> M;
        vector G(N, vector(N, false));
        for (int i = 0; i < M; i++) {
            cin >> a >> b;
            G[a][b] = G[b][a] = true;
        }
        GeneralMatching<mint> gm(G);
        auto res = gm.solve();
        vector<int> ans(N, -1);
        for (auto [x, y] : res) ans[x] = y, ans[y] = x;
        for (int i = 0; i < N; i++) cout << ans[i] << " ";
        cout << endl;
    }
}


int main() {
    //yosupo_general_matching();
    JOI_2026_Matching();
}