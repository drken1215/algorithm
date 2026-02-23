//
// 半環上の行列 (加法・乗法, 行列累乗)
//   SemiRing は「加法」「乗法」が定義されているクラス。コンストラクタで以下の情報を渡す。
//　　・コンストラクタで、ADD (加法), MUL (乗法), ADD_IDENTITY (加法の単位元), MUL_IDENTITY (乗法の単位元)
//　　・何もしなければ、通常の演算子「+」「*」が呼び出される
//
// verified:
//   AtCoder ABC 445 F - Exactly K Steps 2
//     https://atcoder.jp/contests/abc445/tasks/abc445_f
//


#include <bits/stdc++.h>
using namespace std;


// general semiring matrix (define ADD, MUL, ADD_IDENTITY, MUL_IDENTITY)
template<class SemiRing> struct SemiRingMatrix {
    using FuncOperator = function<SemiRing(SemiRing, SemiRing)>;

    // inner value
    int H, W;
    vector<vector<SemiRing>> val;

    // operators
    SemiRing ADD_IDENTITY = SemiRing(), MUL_IDENTITY = SemiRing(1);
    FuncOperator ADD = [](const SemiRing &a, const SemiRing &b) -> SemiRing { 
        return a + b;
    };
    FuncOperator MUL = [](const SemiRing &a, const SemiRing &b) -> SemiRing { 
        return a * b;
    };
    
    // constructors
    SemiRingMatrix() : H(0), W(0) {}
    SemiRingMatrix(int H, int W) : H(H), W(W), val(H, vector<SemiRing>(W, ADD_IDENTITY)) {}
    SemiRingMatrix(int H, int W, SemiRing v) : H(H), W(W), val(H, vector<SemiRing>(W, v)) {}
    SemiRingMatrix(const SemiRingMatrix &mat) : H(mat.H), W(mat.W), val(mat.val)
        , ADD(mat.ADD), MUL(mat.MUL)
        , ADD_IDENTITY(mat.ADD_IDENTITY), MUL_IDENTITY(mat.MUL_IDENTITY) {}
    SemiRingMatrix(int H, int W, const FuncOperator add, const FuncOperator mul
        , const SemiRing &add_identity, const SemiRing &mul_identity) {
        init(H, W, add, mul, add_identity, mul_identity);
    }
    void init(int h, int w, const SemiRing &x) {
        H = h, W = w;
        val.assign(h, vector<SemiRing>(w, x));
    }
    void init(int h, int w, const FuncOperator add, const FuncOperator mul
    , const SemiRing &add_identity, const SemiRing &mul_identity) {
        H = h, W = w;
        ADD = add, MUL = mul;
        ADD_IDENTITY = add_identity, MUL_IDENTITY = mul_identity;
        val.assign(h, vector<SemiRing>(w, ADD_IDENTITY));
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
    vector<SemiRing>& operator [] (int i) { return val[i]; }
    constexpr vector<SemiRing>& operator [] (int i) const { return val[i]; }
    friend constexpr ostream& operator << (ostream &os, const SemiRingMatrix<SemiRing> &mat) {
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
    constexpr bool operator == (const SemiRingMatrix &r) const {
        return this->val == r.val;
    }
    constexpr bool operator != (const SemiRingMatrix &r) const {
        return this->val != r.val;
    }
    
    // arithmetic operators
    constexpr SemiRingMatrix& operator += (const SemiRingMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        assert(ADD_IDENTITY == r.ADD_IDENTITY), assert(MUL_IDENTITY == r.MUL_IDENTITY);
        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                val[i][j] = ADD(val[i][j], r.val[i][j]);
            }
        }
        return *this;
    }
    constexpr SemiRingMatrix& operator *= (const SemiRing &v) {
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                val[i][j] = MUL(val[i][j], v);
        return *this;
    }
    constexpr SemiRingMatrix& operator *= (const SemiRingMatrix &r) {
        assert(width() == r.height());
        assert(ADD_IDENTITY == r.ADD_IDENTITY), assert(MUL_IDENTITY == r.MUL_IDENTITY);
        SemiRingMatrix<SemiRing> res(height(), r.width(), ADD, MUL, ADD_IDENTITY, MUL_IDENTITY);
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                for (int k = 0; k < width(); ++k)
                    res[i][j] = ADD(res[i][j], MUL(val[i][k], r.val[k][j]));
        return (*this) = res;
    }
    constexpr SemiRingMatrix operator + () const { 
        return SemiRingMatrix(*this);
    }
    constexpr SemiRingMatrix operator + (const SemiRingMatrix &r) const { 
        return SemiRingMatrix(*this) += r;
    }
    constexpr SemiRingMatrix operator * (const SemiRing &v) const { 
        return SemiRingMatrix(*this) *= v;
    }
    constexpr SemiRingMatrix operator * (const SemiRingMatrix &r) const { 
        return SemiRingMatrix(*this) *= r;
    }
    constexpr vector<SemiRing> operator * (const vector<SemiRing> &v) const {
        assert(width() == v.size());
        vector<SemiRing> res(height(), ADD_IDENTITY);
        for (int i = 0; i < height(); i++)
            for (int j = 0; j < width(); j++)
                res[i] = ADD(res[i], MUL(val[i][j], v[j]));
        return res;
    }

    // transpose
    constexpr SemiRingMatrix trans() const {
        SemiRingMatrix<SemiRing> res(height(), width(), ADD, MUL, ADD_IDENTITY, MUL_IDENTITY);
        for (int row = 0; row < width(); row++)
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend constexpr SemiRingMatrix<SemiRing> trans(const SemiRingMatrix<SemiRing> &mat) {
        return mat.trans();
    }
    
    // pow
    constexpr SemiRingMatrix pow(long long n) const {
        assert(height() == width());
        SemiRingMatrix<SemiRing> res(height(), width(), ADD, MUL, ADD_IDENTITY, MUL_IDENTITY);
        SemiRingMatrix<SemiRing> mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = MUL_IDENTITY;
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend constexpr SemiRingMatrix<SemiRing> pow(const SemiRingMatrix<SemiRing> &mat, long long n) {
        return mat.pow(n);
    }
};


//------------------------------//
// Examples
//------------------------------//

// AtCoder ABC 445 F - Exactly K Steps 2
void ABC_445_F() {
    long long N, K;
    cin >> N >> K;
    const long long INF = 1LL << 60;
    auto add = [&](long long a, long long b) -> long long { return min(a, b); };
    auto mul = [&](long long a, long long b) -> long long { return a + b; };
    SemiRingMatrix<long long> C(N, N, add, mul, INF, 0);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) cin >> C[j][i];
    auto P = pow(C, K);    
    for (int s = 0; s < N; s++) cout << P[s][s] << endl;
}


int main() {
    ABC_445_F();
}

