//
// mod. p 行列 (行列累乗と、掃き出し法)
//
// verified:
//   TCO 2013 Round 2A Med TheMagicMatrix
//     https://community.topcoder.com/stat?c=problem_statement&pm=12495&rd=15594
// 


#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
using namespace std;


const long long MOD = 1234567891LL;
long long modinv(long long a, long long mod) {
    long long b = mod, u = 1, v = 0;
    while (b) {
        long long t = a/b;
        a -= t*b; swap(a, b);
        u -= t*v; swap(u, v);
    }
    u %= mod;
    if (u < 0) u += mod;
    return u;
}

// matrix
template<int MOD> struct Matrix {
    vector<vector<long long> > val;
    Matrix(int n, int m, long long x = 0) : val(n, vector<long long>(m, x)) {}
    void init(int n, int m, long long x = 0) {val.assign(n, vector<long long>(m, x));}
    size_t size() const {return val.size();}
    inline vector<long long>& operator [] (int i) {return val[i];}
};

template<int MOD> ostream& operator << (ostream& s, Matrix<MOD> A) {
    s << endl; 
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            s << A[i][j] << ", ";
        }
        s << endl;
    }
    return s;
}

template<int MOD> Matrix<MOD> operator * (Matrix<MOD> A, Matrix<MOD> B) {
    Matrix<MOD> R(A.size(), B[0].size());
    for (int i = 0; i < A.size(); ++i) 
        for (int j = 0; j < B[0].size(); ++j)
            for (int k = 0; k < B.size(); ++k) 
                R[i][j] = (R[i][j] + A[i][k] * B[k][j] % MOD) % MOD; 
    return R;
}

template<int MOD> Matrix<MOD> pow(Matrix<MOD> A, long long n) {
    Matrix<MOD> R(A.size(), A.size());
    for (int i = 0; i < A.size(); ++i) R[i][i] = 1;
    while (n > 0) {
        if (n & 1) R = R * A;
        A = A * A;
        n >>= 1;
    }
    return R;
}

template<int MOD> int GaussJordan(Matrix<MOD> &A, bool is_extended = false) {
    int m = A.size(), n = A[0].size();
    for (int row = 0; row < m; ++row)
        for (int col = 0; col < n; ++col)
            A[row][col] = (A[row][col] % MOD + MOD) % MOD;
    int rank = 0;
	for (int col = 0; col < n; ++col) {
        if (is_extended && col == n-1) break;
		int pivot = -1;
        for (int row = rank; row < m; ++row) {
            if (A[row][col] != 0) {
                pivot = row;
                break;
            }
        }
		if (pivot == -1) continue;
        swap(A[pivot], A[rank]);
        auto inv = modinv(A[rank][col], MOD);
        for (int col2 = 0; col2 < n; ++col2)
            A[rank][col2] = A[rank][col2] * inv % MOD;
        for (int row = 0; row < m; ++row) {
            if (row != rank && A[row][col]) {
                auto fac = A[row][col];
                for (int col2 = 0; col2 < n; ++col2) {
                    A[row][col2] -= A[rank][col2] * fac % MOD;
                    if (A[row][col2] < 0) A[row][col2] += MOD;
                }
            }
        }
        ++rank;
    }
    return rank;
}

template<int MOD> int linear_equation(Matrix<MOD> A, vector<long long> b, vector<long long> &res) {
    int m = A.size(), n = A[0].size();
    Matrix<MOD> M(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }
    int rank = GaussJordan(M, true);

    // check if it has no solution
    for (int row = rank; row < m; ++row) if (M[row][n]) return -1;

    // answer
    res.assign(n, 0);
    for (int i = 0; i < rank; ++i) res[i] = M[i][n];
    return rank;
}



long long modpow(long long a, long long n, long long mod) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}

class TheMagicMatrix {
public:
    int find(int n, vector <int> rows, vector <int> cols, vector <int> vals) {
        int m = rows.size();

        // 数字のない行と列がある場合
        set<int> sr, sc;
        for (int i = 0; i < rows.size(); ++i) {
            sr.insert(rows[i]);
            sc.insert(cols[i]);
        }
        if (sr.size() < n && sc.size() < n) {
            long long ex = (n-1)*(n-1) - m + 1;
            return modpow(10LL, ex, MOD);
        }

        // 連立一次方程式を立てる
        Matrix<2> M2(n*2+m, n*n);
        Matrix<5> M5(n*2+m, n*n);
        vector<long long> b(n*2+m);

        // 行和
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                int id = i * n + j;
                M2[i][id] = M5[i][id] = 1;
            }
        }

        // 列和
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                int id = i * n + j;
                M2[j + n][id] = M5[j + n][id] = 1;
            }
        }

        // 条件
        for (int k = 0; k < m; ++k) {
            int id = rows[k] * n + cols[k];
            M2[k + n*2][id] = M5[k + n*2][id] = 1;
            b[k + n*2] = vals[k];
        }

        // X = 0, 1, ..., 9
        long long res = 0;
        for (int X = 0; X < 10; ++X) {
            for (int i = 0; i < n*2; ++i) b[i] = X;

            vector<long long> v;
            int rank2 = linear_equation(M2, b, v);
            int rank5 = linear_equation(M5, b, v);
            if (rank2 == -1 || rank5 == -1) continue;
            long long tmp = modpow(2LL, n*n-rank2, MOD) * modpow(5LL, n*n-rank5, MOD) % MOD;
            res = (res + tmp) % MOD;
        }
        return res;
    }
};






int main() {
    int n, m;
    while (cin >> n >> m) {
        vector<int> r(m), c(m), v(m);
        for (int i = 0; i < m; ++i) cin >> r[i] >> c[i] >> v[i];
        TheMagicMatrix tmm;
        cout << tmm.find(n, r, c, v) << endl;
    }
}
