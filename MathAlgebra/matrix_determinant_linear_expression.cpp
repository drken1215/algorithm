//
// mod. p 行列 (行列累乗、掃き出し法)
//
// verified:
//   AtCoder ARC 176 D - Swap Permutation
//     https://atcoder.jp/contests/arc176/tasks/arc176_d
//
//   TCO 2013 Round 2A Med TheMagicMatrix
//     https://vjudge.net/problem/TopCoder-12495
//
//   AOJ 3369 (?) Namori Counting (OUPC 2023 day2-D)
//     https://onlinejudge.u-aizu.ac.jp/beta/room.html#OUPC2023Day2/problems/D
//
//   ARC 199 B - Adjacent Replace
//     https://atcoder.jp/contests/arc199/tasks/arc199_b
//


#include <bits/stdc++.h>
using namespace std;


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



//------------------------------//
// Examples
//------------------------------//

// AtCoder ARC 176 D - Swap Permutation
void ARC_176_D() {
    const int MOD = 998244353;
    using mint = Fp<MOD>;
    
    long long N, M;
    cin >> N >> M;
    vector<long long> P(N);
    for (int i = 0; i < N; ++i) cin >> P[i], --P[i];
    
    long long NC = N * (N - 1) / 2, N2C = (N - 2) * (N - 3) / 2;
    mint all = mint(NC).pow(M);
    
    if (N == 2) {
        cout << all << endl;
        return;
    }
    
    auto Q = P;
    sort(Q.begin(), Q.end());
    vector<long long> left(N+1, 0), right(N+1, 0);
    for (int i = 0; i < N; ++i) {
        left[i+1] = left[i] + Q[i];
        right[i+1] = right[i] + Q[N-i-1];
    }
    
    auto calc_sum = [&](long long x) -> long long {
        long long l = lower_bound(Q.begin(), Q.end(), x) - Q.begin();
        return (x * l - left[l]) + (right[N - l] - x * (N - l));
    };
    
    vector<mint> f(N, 0);
    mint S = 0;
    for (int i = 0; i < N; ++i) {
        f[i] = calc_sum(P[i]);
        S += f[i];
    }
    
    MintMatrix<mint> A(4, 4);
    A[0][0] = mint(N2C + 1); // / NC;
    A[0][1] = A[0][2] = A[1][2] = A[2][1] = mint(1); // / NC;
    A[1][0] = A[2][0] = mint(N - 2); // / NC;
    A[1][1] = A[2][2] = mint(N2C + N - 2); // / NC;
    A[1][3] = A[2][3] = mint(2); // / NC;
    A[3][1] = A[3][2] = mint(N - 3); // / NC;
    A[3][3] = mint(N2C + N * 2 - 7); // / NC;
    auto AM = pow(A, M);
    
    mint res = 0;
    for (int i = 0; i + 1 < N; ++i) {
        mint diff = abs(P[i] - P[i+1]);
        res += AM[0][0] * diff;
        res += AM[1][0] * (f[i] - diff) / mint(N - 2);
        res += AM[2][0] * (f[i+1] - diff) / mint(N - 2);
        if (N > 3) res += AM[3][0] * (mint(S) / 2 - f[i] - f[i+1] + diff) / mint(N2C);
    }
    cout << res << endl;
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
    //ARC_176_D();
    //AOJ_3369();
    //TCO_2013_Round2_A();
    ARC_199_B();
}