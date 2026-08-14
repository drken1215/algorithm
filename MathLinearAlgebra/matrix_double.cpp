//
// 実数行列 (行列累乗と、掃き出し法)
//
// verified:
//   AOJ 2171 Strange Couple (for gauss-jordan)
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2171
//
//   AOJ 1328 Find the Outlier (for gauss-jordan)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1328
//
//   数学アルゴ本 100 - Simulation of Chemicals (for pow)
//     https://atcoder.jp/contests/math-and-algorithm/tasks/math_and_algorithm_bv
//


#include <bits/stdc++.h>
using namespace std;


// double matrix
template<class DD> struct DoubleMatrix {
    // basic settings
    DD EPS = 1e-10;  // to be set appropriately

    // inner value
    int H, W;
    vector<vector<DD>> val;
    
    // constructors
    DoubleMatrix() {}
    DoubleMatrix(const DoubleMatrix&) = default;
    DoubleMatrix& operator = (const DoubleMatrix&) = default;
    DoubleMatrix(int h, int w) : H(h), W(w), val(h, vector<DD>(w)) {}
    DoubleMatrix(int h, int w, DD x) : H(h), W(w), val(h, vector<DD>(w, x)) {}
    void init(int h, int w, DD x) {
        H = h, W = w;
        val.assign(h, vector<DD>(w, x));
    }
    void resize(int h, int w) {
        H = h, W = w;
        val.resize(h);
        for (int i = 0; i < h; ++i) val[i].resize(w);
    }
    constexpr DD get_eps() { return EPS; }
    constexpr void set_eps(DD eps) { EPS = eps; }
    
    // getter and debugger
    constexpr int height() const { return H; }
    constexpr int width() const { return W; }
    constexpr bool empty() const { return height() == 0; }
    vector<DD>& operator [] (int i) { return val[i]; }
    const vector<DD>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const DoubleMatrix &mat) {
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
    constexpr bool operator == (const DoubleMatrix &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const DoubleMatrix &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr DoubleMatrix& operator += (const DoubleMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] + r.val[i][j];
        return *this;
    }
    constexpr DoubleMatrix& operator -= (const DoubleMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] - r.val[i][j];
        return *this;
    }
    constexpr DoubleMatrix& operator *= (const DD &v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = val[i][j] * v;
        return *this;
    }
    constexpr DoubleMatrix& operator *= (const DoubleMatrix &r) {
        assert(width() == r.height());
        DoubleMatrix<DD> res(height(), r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] = res[i][j] + val[i][k] * r.val[k][j];
        return (*this) = res;
    }
    constexpr DoubleMatrix operator + () const { 
        return DoubleMatrix(*this);
    }
    constexpr DoubleMatrix operator - () const {
        DoubleMatrix res(*this);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res.val[i][j] = -res.val[i][j];
        return res;
    }
    constexpr DoubleMatrix operator + (const DoubleMatrix &r) const { 
        return DoubleMatrix(*this) += r;
    }
    constexpr DoubleMatrix operator - (const DoubleMatrix &r) const {
        return DoubleMatrix(*this) -= r;
    }
    constexpr DoubleMatrix operator * (const DD &v) const {
        return DoubleMatrix(*this) *= v;
    }
    constexpr DoubleMatrix operator * (const DoubleMatrix &r) const {
        return DoubleMatrix(*this) *= r;
    }
    constexpr vector<DD> operator * (const vector<DD> &v) const {
        assert(width() == v.size());
        vector<DD> res(height(), DD(0));
        for (int i = 0; i < height(); i++)
            for (int j = 0; j < width(); j++)
                res[i] += val[i][j] * v[j];
        return res;
    }

    // transpose
    constexpr DoubleMatrix trans() const {
        DoubleMatrix<DD> res(width(), height());
        for (int row = 0; row < width(); row++)
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend constexpr DoubleMatrix trans(const DoubleMatrix &mat) {
        return mat.trans();
    }
    
    // pow
    constexpr DoubleMatrix pow(long long n) const {
        assert(height() == width());
        DoubleMatrix<DD> res(height(), width());
        DoubleMatrix<DD> mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = DD(1);
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr DoubleMatrix pow(const DoubleMatrix &mat, long long n) {
        return mat.pow(n);
    }
    
    // gauss-jordan
    constexpr int find_pivot(int cur_rank, int col) const {
        int pivot = -1;
        DD max_v = EPS;
        for (int row = cur_rank; row < height(); ++row) {
            if (abs(val[row][col]) > max_v) {
                max_v = abs(val[row][col]);
                pivot = row;
            }
        }
        return pivot;
    }
    constexpr void sweep(int cur_rank, int col, int pivot) {
        swap(val[pivot], val[cur_rank]);
        auto fac = val[cur_rank][col];
        for (int col2 = 0; col2 < width(); ++col2) {
            val[cur_rank][col2] /= fac;
        }
        for (int row = 0; row < height(); ++row) {
            if (row != cur_rank && abs(val[row][col]) > EPS) {
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
    friend constexpr int gauss_jordan(DoubleMatrix<DD> &mat, int not_sweep_width = 0) {
        return mat.gauss_jordan(not_sweep_width);
    }

    // rank
    constexpr int get_rank() const {
        if (height() == 0 || width() == 0) return 0;
        DoubleMatrix A(*this);
        if (height() < width()) A = A.trans();
        return A.gauss_jordan(0, false);
    }
    friend constexpr int get_rank(const DoubleMatrix &mat) {
        return mat.get_rank();
    }

    // find one solution
    friend constexpr int linear_equation
    (const DoubleMatrix &mat, const vector<DD> &b, vector<DD> &res) {
        // extend
        DoubleMatrix<DD> A(mat.height(), mat.width() + 1);
        A.set_eps(mat.EPS);
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) A[i][j] = mat.val[i][j];
            A[i].back() = b[i];
        }
        int rank = A.gauss_jordan(1);
        
        // check if it has no solution
        for (int row = rank; row < mat.height(); ++row) if (abs(A[row].back()) > A.EPS) return -1;

        // answer
        res.assign(mat.width(), 0);
        for (int i = 0; i < rank; ++i) res[i] = A[i].back();
        return rank;
    }
    friend constexpr int linear_equation(const DoubleMatrix &mat, const vector<DD> &b) {
        vector<DD> res;
        return linear_equation(mat, b, res);
    }

    // find all solutions
    friend int linear_equation
    (const DoubleMatrix &mat, const vector<DD> &b, vector<DD> &res, vector<vector<DD>> &zeros) {
        // extend
        DoubleMatrix<DD> A(mat.height(), mat.width() + 1);
        A.set_eps(mat.EPS);
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) A[i][j] = mat.val[i][j];
            A[i].back() = b[i];
        }
        vector<int> core;
        int rank = A.gauss_jordan(core, 1);
        
        // check if it has no solution
        for (int row = rank; row < mat.height(); ++row) {
            if (abs(A[row].back()) > A.EPS) return -1;
        }

        // construct the core solution
        res.assign(mat.width(), DD(0));
        for (int i = 0; i < (int)core.size(); i++) res[core[i]] = A[i].back();
    
        // construct the all solutions
        zeros.clear();
        vector<bool> use(mat.width(), 0);
        for (auto c : core) use[c] = true;
        for (int j = 0; j < mat.width(); j++) {
            if (use[j]) continue;
            vector<DD> zero(mat.width(), DD(0));
            zero[j] = DD(1);
            for (int i = 0; i < (int)core.size(); i++) zero[core[i]] = -A[i][j];
            zeros.push_back(zero);
        }
        return rank;
    }

    // determinant
    constexpr DD det() const {
        assert(height() == width());
        if (height() == 0) return DD(1);
        DoubleMatrix<DD> A(*this);
        int rank = 0;
        DD res = DD(1);
        for (int col = 0; col < width(); ++col) {
            int pivot = A.find_pivot(rank, col);
            if (pivot == -1) return DD(0);
            if (pivot != rank) res = -res;
            res *= A[pivot][rank];
            A.sweep(rank++, col, pivot, false);
        }
        return res;
    }
    friend constexpr DD det(const DoubleMatrix &mat) {
        return mat.det();
    }

    // inv
    constexpr DoubleMatrix inv() const {
        assert(height() == width());

        // extend
        DoubleMatrix<DD> A(height(), width() + height());
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) A[i][j] = val[i][j];
            A[i][i+width()] = DD(1);
        }
        vector<int> core;
        int rank = A.gauss_jordan(height(), true);

        // gauss jordan
        if (rank < height()) return DoubleMatrix();
        DoubleMatrix<DD> res(height(), width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                res[i][j] = A[i][j+width()];
        return res;
    }
    friend constexpr DoubleMatrix inv(const DoubleMatrix &mat) {
        return mat.inv();
    }
};



//------------------------------//
// Examples
//------------------------------//

// AOJ 2171 Strange Couple
void AOJ_2171() {
    int N, s, t;
    while (cin >> N >> s >> t, N) {
        --s, --t;
        vector<int> q(N);
        vector<vector<int>> a(N, vector<int>(N));
        for (int i = 0; i < N; ++i) cin >> q[i];
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) cin >> a[i][j];

        // Dijkstra
        const int INF = 1<<29;
        vector<int> dist(N, INF);
        vector<bool> seen(N, 0);
        dist[t] = 0;
        for (int iter = 0; iter < N; ++iter) {
            int curd = INF;
            int v = -1;
            for (int i = 0; i < N; ++i) {
                if (seen[i]) continue;
                if (curd > dist[i]) {
                    curd = dist[i];
                    v = i;
                }
            }
            if (v == -1) break;
            for (int w = 0; w < N; ++w) {
                if (w == v) continue;
                if (a[v][w] == 0) continue;
                dist[w] = min(dist[w], curd + a[v][w]);
            }
            seen[v] = true;
        }
        if (dist[s] >= INF) {
            cout << "impossible" << endl;
            return;
        }

        // 連立一次方程式を作る
        DoubleMatrix<double> A(N, N, 0);
        vector<double> b(N, 0);
        for (int v = 0; v < N; ++v) {
            if (v == t) {
                A[v][v] = 1;
                b[v] = 0;
            }
            else {
                vector<int> neigbor;
                for (int w = 0; w < N; ++w) {
                    if (a[v][w] == 0) continue;
                    if (q[v] == 1 && dist[w] + a[v][w] != dist[v]) continue;
                    neigbor.push_back(w);
                }
                int K = neigbor.size();
                for (auto w : neigbor) {
                    A[v][w] -= 1;
                    b[v] += a[v][w];
                }
                A[v][v] += K;
            }
        }
          
        // 解く
        vector<double> res;
        auto rank = linear_equation(A, b, res);
        if (res.empty()) cout << "impossible" << endl;
        else cout << fixed << setprecision(15) << res[s] << endl;
    }
}

// AOJ 1328 Find the Outlier
void AOJ_1328() {
    using D = long double;
    
    auto dpow = [&](D a, int n) -> D {
        D res = 1.0;
        for (int i = 0; i < n; ++i) res *= a;
        return res;
    };
    auto func = [&](const vector<D> &coef, int i) -> D {
        D res = 0.0;
        for (int p = 0; p < (int)coef.size(); ++p)
            res += coef[p] * pow(i, p);
        return res;
    };
    
    int d;
    while (cin >> d, d) {
        vector<D> v(d + 3);
        for (int i = 0; i < d + 3; ++i) cin >> v[i];

        bool finish = false;
        int res = 0;
        for (int i = 0; i < d + 3 && !finish; ++i) {
            for (int j = i + 1; j < d + 3 && !finish; ++j) {
                DoubleMatrix<D> A(d + 1, d + 1);
                A.set_eps(1e-5);
                vector<D> b(d + 1);
                for (int k = 0, iter = 0; k < d+3; ++k) {
                    if (k == i || k == j) continue;
                    for (int p = 0; p < d + 1; ++p) {
                        A[iter][p] = pow(k, p);
                        b[iter] = v[k];
                    }
                    ++iter;
                }
                vector<D> ans;
                auto rank = linear_equation(A, b, ans);
                if (ans.empty()) continue;
                D vi = func(ans, i), vj = func(ans, j);
                int num = 0;
                if (fabs(vi - v[i]) > A.EPS) res = i, ++num;
                if (fabs(vj - v[j]) > A.EPS) res = j, ++num;
                if (num == 1) goto end;
            }
        }
        end:
        cout << res << endl;
    }
}

// 数学アルゴ本 100 - Simulation of Chemicals
void MathAlgo100() {
    using DD = long double;
    int Q;
    cin >> Q;
    while (Q--) {
        DD X, Y, Z;
        long long T;
        cin >> X >> Y >> Z >> T;

        DoubleMatrix<DD> A(3, 3, 0);
        A[0][0] = 1.0 - X, A[0][1] = Y;
        A[1][1] = 1.0 - Y, A[1][2] = Z;
        A[2][2] = 1.0 - Z, A[2][0] = X;
        auto P = A.pow(T);
        vector<DD> ini({1, 1, 1});
        auto res = P * ini;
        cout << fixed << setprecision(10);
        cout << res[0] << " " << res[1] << " " << res[2] << endl;
    }
}


int main() {
    AOJ_2171();
    //AOJ_1328();
    //MathAlgo100();
}