//
// 実数行列 (行列累乗と、掃き出し法)
//
// verified:
//   AOJ 2171 Strange Couple
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2171
// 


#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl


using D = double;
const D EPS = 1e-10;

template<class T> struct Matrix {
    vector<vector<T> > val;
    Matrix(int n, int m, T x = 0) : val(n, vector<T>(m, x)) {}
    void init(int n, int m, T x = 0) {val.assign(n, vector<T>(m, x));}
    size_t size() const {return val.size();}
    inline vector<T>& operator [] (int i) {return val[i];}
};

template<class T> ostream& operator << (ostream& s, Matrix<T> A) {
    s << endl; 
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            s << A[i][j] << ", ";
        }
        s << endl;
    }
    return s;
}

template<class T> Matrix<T> operator * (Matrix<T> A, Matrix<T> B) {
    Matrix<T> R(A.size(), B[0].size());
    for (int i = 0; i < A.size(); ++i)
        for (int j = 0; j < B[0].size(); ++j)
            for (int k = 0; k < B.size(); ++k)
                R[i][j] += A[i][k] * B[k][j];
    return R;
}

template<class T> Matrix<T> pow(Matrix<T> A, long long n) {
    Matrix<T> R(A.size(), A.size());
    for (int i = 0; i < A.size(); ++i) R[i][i] = 1;
    while (n > 0) {
        if (n & 1) R = R * A;
        A = A * A;
        n >>= 1;
    }
    return R;
}

template<class T> int GaussJordan(Matrix<T> &A) {
    int rank = 0;
    for (int col = 0; col < A[0].size(); ++col) {
        int pivot = -1;
        T ma = EPS;
        for (int row = rank; row < A.size(); ++row) {
            if (abs(A[row][col]) > ma) {
                ma = abs(A[row][col]);
                pivot = row;
            }
        }
        if (pivot == -1) continue;
        swap(A[pivot], A[rank]);
        auto fac = A[rank][col];
        for (int col2 = 0; col2 < A[0].size(); ++col2) A[rank][col2] /= fac;
        for (int row = 0; row < A.size(); ++row) {
            if (row != rank && abs(A[row][col]) > EPS) {
                auto fac = A[row][col];
                for (int col2 = 0; col2 < A[0].size(); ++col2) {
                    A[row][col2] -= A[rank][col2] * fac;
                }
            }
        }
        ++rank;
    }
    return rank;
}

template<class T> vector<T> linear_equation(Matrix<T> A, vector<T> b) {
    Matrix<T> M(A.size(), A[0].size() + 1);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[0].size(); ++j) M[i][j] = A[i][j];
        M[i][(int)A[0].size()] = b[i];
    }
    int rank = GaussJordan(M);

    // check if it has no solution
    vector<T> res;
    if (rank > 0) {
        bool exist = false;
        for (int col = 0; col < A[0].size(); ++col)
            if (abs(M[rank-1][col]) > EPS) exist = true;
        if (!exist) return res;
    }
    for (int row = rank; row < A.size(); ++row)
        if (abs(M[row][A[0].size()]) > EPS)
            return res;

    // answer
    res.assign(A[0].size(), 0);
    for (int i = 0; i < rank; ++i) res[i] = M[i][A[0].size()];
    return res;
}



int N, s, t;
vector<int> q;
vector<vector<int> > a;

void solve() {
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

    // 連立一次方程式
    Matrix<D> A(N, N, 0); vector<D> b(N, 0);
    for (int v = 0; v < N; ++v) {
        if (v == t) {
            A[v][v] = 1;
            b[v] = 0;
        }
        else if (q[v] == 1) {
            int K = 0;
            for (int w = 0; w < N; ++w) {
                if (a[v][w] == 0) continue;
                if (dist[w] + a[v][w] != dist[v]) continue;
                ++K;
                A[v][w] -= 1;
                b[v] += a[v][w];
            }
            A[v][v] += K;
        }
        else {
            int K = 0;
            for (int w = 0; w < N; ++w) {
                if (a[v][w] == 0) continue;
                ++K;
                A[v][w] -= 1;
                b[v] += a[v][w];
            }
            A[v][v] += K;
        }
    }
    auto res = linear_equation(A, b); 
    if (res.empty()) cout << "impossible" << endl;
    else cout << fixed << setprecision(15) << res[s] << endl;
}

int main() {
    while (cin >> N >> s >> t, N) {
        --s, --t;
        q.resize(N);
        for (int i = 0; i < N; ++i) cin >> q[i];
        a.assign(N, vector<int>(N));
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) cin >> a[i][j];
        solve();
    }
}
