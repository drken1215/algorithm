//
// 行列 (基本演算)
//
// verified:
//   AGC 003 F - Fraction of Fractal
//     https://atcoder.jp/contests/agc003/tasks/agc003_f
//


#include <iostream>
#include <vector>
using namespace std;


// matrix
template<class T> struct Matrix {
    vector<vector<T> > val;
    Matrix(int n = 1, int m = 1, T v = 0) : val(n, vector<T>(m, v)) {}
    void init(int n, int m, T v = 0) {val.assign(n, vector<T>(m, v));}
    void resize(int n, int m) {
        val.resize(n);
        for (int i = 0; i < n; ++i) val[i].resize(m);
    }
    Matrix<T>& operator = (const Matrix<T> &A) {
        val = A.val;
        return *this;
    }
    size_t size() const {return val.size();}
    vector<T>& operator [] (int i) {return val[i];}
    const vector<T>& operator [] (int i) const {return val[i];}
    friend ostream& operator << (ostream& s, const Matrix<T>& M) {
        s << endl;
        for (int i = 0; i < (int)M.size(); ++i) s << M[i] << endl;
        return s;
    }
};

template<class T> Matrix<T> operator * (const Matrix<T> &A, const Matrix<T> &B) {
    Matrix<T> R(A.size(), B[0].size());
    for (int i = 0; i < A.size(); ++i)
        for (int j = 0; j < B[0].size(); ++j)
            for (int k = 0; k < B.size(); ++k)
                R[i][j] += A[i][k] * B[k][j];
    return R;
}

template<class T> Matrix<T> pow(const Matrix<T> &A, long long n) {
    Matrix<T> R(A.size(), A.size());
    auto B = A;
    for (int i = 0; i < A.size(); ++i) R[i][i] = 1;
    while (n > 0) {
        if (n & 1) R = R * B;
        B = B * B;
        n >>= 1;
    }
    return R;
}

template<class T> vector<T> operator * (const Matrix<T> &A, const vector<T> &B) {
    vector<T> v(A.size());
    for (int i = 0; i < A.size(); ++i)
        for (int k = 0; k < B.size(); ++k)
            v[i] += A[i][k] * B[k];
    return v;
}

template<class T> Matrix<T> operator + (const Matrix<T> &A, const Matrix<T> &B) {
    Matrix<T> R(A.size(), A[0].size());
    for (int i = 0; i < A.size(); ++i)
        for (int j = 0; j < A[0].size(); ++j)
            R[i][j] = A[i][j] + B[i][j];
    return R;
}

template<class T> Matrix<T> operator - (const Matrix<T> &A, const Matrix<T> &B) {
    Matrix<T> R(A.size(), A[0].size());
    for (int i = 0; i < A.size(); ++i)
        for (int j = 0; j < A[0].size(); ++j)
            R[i][j] = A[i][j] - B[i][j];
    return R;
}


// modint
template<int MOD> struct Fp {
    long long val;
    constexpr Fp(long long v = 0) noexcept : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr int getmod() { return MOD; }
    constexpr Fp operator - () const noexcept {
        return val ? MOD - val : 0;
    }
    constexpr Fp operator + (const Fp& r) const noexcept { return Fp(*this) += r; }
    constexpr Fp operator - (const Fp& r) const noexcept { return Fp(*this) -= r; }
    constexpr Fp operator * (const Fp& r) const noexcept { return Fp(*this) *= r; }
    constexpr Fp operator / (const Fp& r) const noexcept { return Fp(*this) /= r; }
    constexpr Fp& operator += (const Fp& r) noexcept {
        val += r.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    constexpr Fp& operator -= (const Fp& r) noexcept {
        val -= r.val;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr Fp& operator *= (const Fp& r) noexcept {
        val = val * r.val % MOD;
        return *this;
    }
    constexpr Fp& operator /= (const Fp& r) noexcept {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b; swap(a, b);
            u -= t * v; swap(u, v);
        }
        val = val * u % MOD;
        if (val < 0) val += MOD;
        return *this;
    }
    constexpr bool operator == (const Fp& r) const noexcept {
        return this->val == r.val;
    }
    constexpr bool operator != (const Fp& r) const noexcept {
        return this->val != r.val;
    }
    friend constexpr ostream& operator << (ostream &os, const Fp<MOD>& x) noexcept {
        return os << x.val;
    }
    friend constexpr Fp<MOD> modpow(const Fp<MOD> &a, long long n) noexcept {
        if (n == 0) return 1;
        auto t = modpow(a, n / 2);
        t = t * t;
        if (n & 1) t = t * a;
        return t;
    }
};


const int MOD = 1000000007;
using mint = Fp<MOD>;
long long H, W, K;
vector<string> fi;

mint solve() {
    if (K <= 1) return 1;
    long long N = 0, yoko = 0, tate = 0, adj = 0;
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j)
            if (fi[i][j] == '#') ++N;
    for (int i = 0; i < H; ++i)
        if (fi[i][0] == '#' && fi[i].back() == '#') ++yoko;
    for (int j = 0; j < W; ++j)
        if (fi[0][j] == '#' && fi[H-1][j] == '#') ++tate;
    if (yoko && tate) return 1;
    if (!yoko && !tate) return modpow(mint(N), K-1);

    Matrix<mint> A(2, 2, 0);
    A[0][0] = N;
    if (yoko) {
        for (int i = 0; i < H; ++i)
            for (int j = 0; j + 1 < W; ++j)
                if (fi[i][j] == '#' && fi[i][j+1] == '#') ++adj;
        A[0][1] = -adj;
        A[1][1] = yoko;
    }
    else {
        for (int i = 0; i + 1 < H; ++i)
            for (int j = 0; j < W; ++j)
                if (fi[i][j] == '#' && fi[i+1][j] == '#') ++adj;
        A[0][1] = -adj;
        A[1][1] = tate;
    }
    auto res = pow(A, K-1);
    return res[0][0] + res[0][1];
}

int main() {
    cin >> H >> W >> K;
    fi.resize(H);
    for (int i = 0; i < H; ++i) cin >> fi[i];
    cout << solve() << endl;
}
