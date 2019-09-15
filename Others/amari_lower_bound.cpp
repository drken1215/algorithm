//
// pn + r (n は非負整数) で表せる整数のうち、x 以上となる最小の整数
//
// verified:
//   ABC 129 F - Takahashi's Basics in Education and Learning
//     https://atcoder.jp/contests/abc129/tasks/abc129_f
//


#include <iostream>
#include <vector>
using namespace std;



// pn + r (n は非負整数) で表せる整数のうち、x 以上となる最小の整数
long long amari_lower_bound(long long p, long long r, long long x) {
    if (r >= x) return r;
    return (x - r + p-1) / p * p + r;
}



// matrix
int MOD;
struct Matrix {
    vector<vector<long long> > val;
    Matrix(int n, int m, long long x = 0) : val(n, vector<long long>(m, x)) {}
    void init(int n, int m, long long x = 0) {val.assign(n, vector<long long>(m, x));}
    size_t size() const {return val.size();}
    inline vector<long long>& operator [] (int i) {return val[i];}
};

Matrix operator * (Matrix A, Matrix B) {
    Matrix R(A.size(), B[0].size());
    for (int i = 0; i < A.size(); ++i) 
        for (int j = 0; j < B[0].size(); ++j)
            for (int k = 0; k < B.size(); ++k) 
                R[i][j] = (R[i][j] + A[i][k] * B[k][j] % MOD) % MOD; 
    return R;
}

Matrix pow(Matrix A, long long n) {
    Matrix R(A.size(), A.size());
    for (int i = 0; i < A.size(); ++i) R[i][i] = 1;
    while (n > 0) {
        if (n & 1) R = R * A;
        A = A * A;
        n >>= 1;
    }
    return R;
}

long long modpow(long long a, long long n) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % MOD;
        a = a * a % MOD;
        n >>= 1;
    }
    return res;
}

// x^(n-1) + x^(n-2) + ... + 1
long long ser(long long x, long long n) {
    if (n == 0) return 0;
    Matrix M(2, 2);
    M[0][0] = x % MOD; M[0][1] = 1;
    M[1][0] = 0; M[1][1] = 1;
    auto P = pow(M, n-1);
    return (P[0][0] + P[0][1]) % MOD;
}

// x^(n-1) + 2x^(n-2) + 3x^(n-3) + ... + n
long long ser2(long long x, long long n) {
    if (n == 0) return 0;
    Matrix M(3, 3);
    M[0][0] = x % MOD; M[0][1] = 1; M[0][2] = 1;
    M[1][0] = 0; M[1][1] = 1; M[1][2] = 1;
    M[2][0] = 0; M[2][1] = 0; M[2][2] = 1;
    auto P = pow(M, n-1);
    return (P[0][0] + P[0][1] + P[0][2]) % MOD;
}

long long solve(long long L, long long A, long long B, long long M) {
    MOD = M;
    long long C = A + B * (L-1);
    long long res = 0;
    long long curd = 0;
    long long sum = 0;
    long long fac = 1;
    for (long long d = 18; d >= 1; --d) {
        // まずは d 桁の個数を数える
        long long beki = 1;
        for (int i = 0; i < d - 1; ++i) beki *= 10;
        long long low = beki, high = beki * 10 - 1;

        high = min(high, C);
        if (high < A || C < low) continue;

        long long low_kou = amari_lower_bound(B, A, low);
        long long up_kou = amari_lower_bound(B, A, high + 1);
        long long num = (up_kou - low_kou) / B;

        // それを使って計算
        long long x = modpow(10LL, d);
        long long alr = num % MOD * sum % MOD;
        long long add = ser2(x, num) * fac % MOD;
        res = (res + (add + alr) % MOD * (B % MOD) % MOD) % MOD;
        sum = (sum + fac * ser(x, num) % MOD) % MOD;
        fac = fac * modpow(modpow(10LL, d), num) % MOD;
    }
    long long hosei = ((A - B) % M + M) % M * sum % M;
    res = (res + hosei) % M;
    return res;   
}

int main() {
    long long L, A, B, M;
    while (cin >> L >> A >> B >> M) cout << solve(L, A, B, M) << endl;
}
