//
// binary 行列 (行列累乗と、掃き出し法)
//
// verified:
//   みんなのプロコン 2019 E - Odd Subrectangles
//     https://atcoder.jp/contests/yahoo-procon2019-qual/tasks/yahoo_procon2019_qual_e
//
//   AOJ 1308 Awkward Lights (ICPC アジア 2010 D)
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1308
//
//   AOJ 2624 Graph Automata Player
//     judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2624
//


#include <iostream>
#include <vector>
#include <bitset>
using namespace std;


const int MAX_ROW = 510; // to be set appropriately
const int MAX_COL = 510; // to be set appropriately
struct BitMatrix {
    int H, W;
    bitset<MAX_COL> val[MAX_ROW];
    BitMatrix(int m = 1, int n = 1) : H(m), W(n) {}
    inline bitset<MAX_COL>& operator [] (int i) {return val[i];}
};

ostream& operator << (ostream& s, BitMatrix A) {
    s << endl;
    for (int i = 0; i < A.H; ++i) {
        for (int j = 0; j < A.W; ++j) {
            s << A[i][j] << ", ";
        }
        s << endl;
    }
    return s;
}

inline BitMatrix operator * (BitMatrix A, BitMatrix B) {
    BitMatrix R(A.H, B.W);
    BitMatrix tB(B.W, B.H);
    for (int i = 0; i < tB.H; ++i) for (int j = 0; j < tB.W; ++j) tB[i][j] = B[j][i];
    for (int i = 0; i < R.H; ++i) for (int j = 0; j < R.W; ++j) R[i][j] = ((A[i] & tB[j]).count() & 1);
    return R;
}

inline BitMatrix pow(BitMatrix A, long long n) {
    BitMatrix R(A.H, A.H);
    for (int i = 0; i < A.H; ++i) R[i][i] = 1;
    while (n > 0) {
        if (n & 1) R = R * A;
        A = A * A;
        n >>= 1;
    }
    return R;
}

int GaussJordan(BitMatrix &A, bool is_extended = false) {
    int rank = 0;
    for (int col = 0; col < A.W; ++col) {
        if (is_extended && col == A.W - 1) break;
        int pivot = -1;
        for (int row = rank; row < A.H; ++row) {
            if (A[row][col]) {
                pivot = row;
                break;
            }
        }
        if (pivot == -1) continue;
        swap(A[pivot], A[rank]);
        for (int row = 0; row < A.H; ++row) {
            if (row != rank && A[row][col]) A[row] ^= A[rank];
        }
        ++rank;
    }
    return rank;
}

int linear_equation(BitMatrix A, vector<int> b, vector<int> &res) {
    int m = A.H, n = A.W;
    BitMatrix M(m, n + 1);
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


const int MOD = 998244353;
long long modpow(long long a, long long n, long long mod) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}

int main() {
    int N, M; cin >> N >> M;
    BitMatrix A(N, M);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            int a; cin >> a;
            if (a) A[i][j] = 1;
        }
    }
    vector<int> res;
    int r = GaussJordan(A);
    cout << (modpow(2LL, N+M-1, MOD) - modpow(2LL, N+M-r-1, MOD) + MOD) % MOD << endl;
}
