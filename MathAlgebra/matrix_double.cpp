//
// 実数行列 (行列累乗と、掃き出し法)
//
// verified:
//   AOJ 2171 Strange Couple
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2171
//
//   AOJ 1328 Find the Outlier
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1328
//


#include <bits/stdc++.h>
using namespace std;


// basic settings
long double EPS = 1e-10;  // to be set appropriately

// matrix
template<class T> struct Matrix {
    // inner value
    vector<vector<T>> val;
    
    // constructors
    Matrix(int H, int W, T x = 0) : val(H, vector<T>(W, x)) {}
    Matrix(const Matrix &mat) : val(mat.val) {}
    void init(int H, int W, T x = 0) {
        val.assign(H, vector<T>(W, x));
    }
    void resize(int H, int W) {
        val.resize(H);
        for (int i = 0; i < H; ++i) val[i].resize(W);
    }
    
    // getter and debugger
    constexpr int height() const { return (int)val.size(); }
    constexpr int width() const { return (int)val[0].size(); }
    vector<T>& operator [] (int i) { return val[i]; }
    constexpr vector<T>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const Matrix<T> &mat) {
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
    constexpr bool operator == (const Matrix &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const Matrix &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr Matrix& operator += (const Matrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                val[i][j] += r[i][j];
            }
        }
        return *this;
    }
    constexpr Matrix& operator -= (const Matrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                val[i][j] -= r[i][j];
            }
        }
        return *this;
    }
    constexpr Matrix& operator *= (T v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] *= v;
        return *this;
    }
    constexpr Matrix& operator *= (const Matrix &r) {
        assert(width() == r.height());
        Matrix<T> res(height(), r.width());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] += val[i][k] * r[k][j];
        return (*this) = res;
    }
    constexpr Matrix operator + () const { return Matrix(*this); }
    constexpr Matrix operator - () const { return Matrix(*this) *= T(-1); }
    constexpr Matrix operator + (const Matrix &r) const { return Matrix(*this) += r; }
    constexpr Matrix operator - (const Matrix &r) const { return Matrix(*this) -= r; }
    constexpr Matrix operator * (T v) const { return Matrix(*this) *= v; }
    constexpr Matrix operator * (const Matrix &r) const { return Matrix(*this) *= r; }
    
    // pow
    constexpr Matrix pow(long long n) const {
        assert(height() == width());
        Matrix<T> res(height(), width()),  mul(*this);
        while (n > 0) {
            if (n & 1) res *= mul;
            mul *= mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr Matrix<T> pow(const Matrix<T> &mat, long long n) {
        return mat.pow(n);
    }
    
    // gauss-jordan
    constexpr int find_pivot(int cur_rank, int col) const {
        int pivot = -1;
        T max_v = EPS;
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
    friend constexpr int gauss_jordan(Matrix<T> &mat, int not_sweep_width = 0) {
        return mat.gauss_jordan(not_sweep_width);
    }
    friend constexpr vector<T> linear_equation(const Matrix<T> &mat, const vector<T> &b) {
        // extend
        Matrix<T> A(mat.height(), mat.width() + 1);
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) A[i][j] = mat.val[i][j];
            A[i].back() = b[i];
        }
        int rank = A.gauss_jordan(1);
        
        // check if it has no solution
        vector<T> res;
        for (int row = rank; row < mat.height(); ++row)
            if (abs(A[row].back()) > EPS)
                return res;

        // answer
        res.assign(mat.width(), 0);
        for (int i = 0; i < rank; ++i) res[i] = A[i].back();
        return res;
    }
};



//------------------------------//
//  Examples
//------------------------------//

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
        Matrix<double> A(N, N, 0);
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
        auto res = linear_equation(A, b);
        if (res.empty()) cout << "impossible" << endl;
        else cout << fixed << setprecision(15) << res[s] << endl;
    }
}

void AOJ_1328() {
    using D = long double;
    EPS = 1e-5;  // set EPS
    
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
                Matrix<D> A(d + 1, d + 1);
                vector<D> b(d + 1);
                for (int k = 0, iter = 0; k < d+3; ++k) {
                    if (k == i || k == j) continue;
                    for (int p = 0; p < d + 1; ++p) {
                        A[iter][p] = pow(k, p);
                        b[iter] = v[k];
                    }
                    ++iter;
                }
                vector<D> ans = linear_equation(A, b);
                if (ans.empty()) continue;
                D vi = func(ans, i), vj = func(ans, j);
                int num = 0;
                if (fabs(vi - v[i]) > EPS) res = i, ++num;
                if (fabs(vj - v[j]) > EPS) res = j, ++num;
                if (num == 1) goto end;
            }
        }
        end:
        cout << res << endl;
    }
}


int main() {
    //AOJ_2171();
    AOJ_1328();
}

