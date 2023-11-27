//
// Gauss Integer
//
// reference:
//   drken: ボードゲーム「共円」に学ぶ、ガウス整数 x + yi の世界
//     https://qiita.com/drken/items/336ef288b451e86c15cb
//
// verified:
//   Yosupo Library Checker - Gcd of Gaussian Integers
//     https://judge.yosupo.jp/problem/gcd_of_gaussian_integers
//
//   New Year Contest L - 机のしみ
//     https://atcoder.jp/contests/NYC2015/tasks/nyc2015_12
//


#include <bits/stdc++.h>
using namespace std;


// Gauss Integer
template<class T> struct GaussNumber {
    // inner value
    T x, y;
    
    // constructor
    constexpr GaussNumber() : x(0), y(0) { }
    constexpr GaussNumber(T x) : x(x), y(0) { }
    constexpr GaussNumber(T x, T y) : x(x), y(y) { }
    
    // various functions
    constexpr GaussNumber conj() const { return GaussNumber(x, -y); }
    constexpr T dot(const GaussNumber &r) const { return x * r.x + y * r.y; }
    constexpr T cross(const GaussNumber &r) const { return x * r.y - y * r.x; }
    constexpr T norm() const { return dot(*this); }
    
    // comparison operators
    constexpr bool operator == (const GaussNumber &r) const {
        return x == r.x && y == r.y;
    }
    constexpr bool operator != (const GaussNumber &r) const {
        return x != r.x || y != r.y;
    }
    
    // arithmetic operators
    constexpr GaussNumber &operator += (const GaussNumber &r) {
        x += r.x, y += r.y;
        return *this;
    }
    constexpr GaussNumber &operator -= (const GaussNumber &r) {
        x -= r.x, y -= r.y;
        return *this;
    }
    constexpr GaussNumber &operator *= (const GaussNumber &r) {
        T tx = x * r.x - y * r.y;
        T ty = x * r.y + y * r.x;
        x = tx, y = ty;
        return *this;
    }
    constexpr GaussNumber &operator /= (T r) {
        assert(r != 0);
        if (r > 0) {
            x = div_near(x, r), y = div_near(y, r);
        } else {
            x = div_near(-x, -r), y = div_near(-y ,-r);
        }
        return *this;
    }
    constexpr GaussNumber &operator /= (const GaussNumber &r) {
        assert(r != 0);
        T nr = r.norm();
        return (*this) = ((*this) * r.conj()) / nr;
    }
    constexpr GaussNumber &operator %= (const GaussNumber &r) {
        assert(r != 0);
        GaussNumber q = (*this) / r;
        return (*this) -= q * r;
    }
    constexpr GaussNumber operator + () const {
        return GaussNumber(*this);
    }
    constexpr GaussNumber operator - () const {
        return GaussNumber(-x, -y);
    }
    constexpr GaussNumber operator + (const GaussNumber &r) const {
        return GaussNumber(*this) += r;
    }
    constexpr GaussNumber operator - (const GaussNumber &r) const {
        return GaussNumber(*this) -= r;
    }
    constexpr GaussNumber operator * (const GaussNumber &r) const {
        return GaussNumber(*this) *= r;
    }
    constexpr GaussNumber operator / (T r) const {
        return GaussNumber(*this) /= r;
    }
    constexpr GaussNumber operator / (const GaussNumber &r) const {
        return GaussNumber(*this) /= r;
    }
    constexpr GaussNumber operator % (const GaussNumber &r) const {
        return GaussNumber(*this) %= r;
    }
    
    // friend operators and functions
    friend ostream &operator << (ostream &os, const GaussNumber &r) {
        if (r.x == 0 && r.y == 0) return os << "0";
        else if (r.x == 0) return os << r.y << "i";
        else os << r.x;
        if (r.y == 0) return os;
        else if (r.y > 0) return os << " + " << r.y << "i";
        else return os << " - " << (-r.y) << "i";
    }
    friend constexpr GaussNumber conj(const GaussNumber &p) {
        return p.conj();
    }
    friend constexpr T dot(const GaussNumber &p, const GaussNumber &q) {
        return p.dot(q);
    }
    friend constexpr T cross(const GaussNumber &p, const GaussNumber &q) {
        return p.cross(q);
    }
    friend constexpr T norm(const GaussNumber &p) {
        return p.norm();
    }
    friend GaussNumber gcd(const GaussNumber &x, const GaussNumber &y) {
        if (y == 0) return x;
        else return gcd(y, x % y);
    }
    
private:
    constexpr T div_near(T x, T y) const {
        assert(y != 0);
        T half = y / 2;
        if (x > half) return (x - half - T(1)) / y + T(1);
        else if (x < -half) return (x + half + T(1)) / y - T(1);
        else return T(0);
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

//Yosupo Library Checker - Gcd of Gaussian Integers
void Yosupo_Gcd_of_Gaussian_Integers() {
    using gint = GaussNumber<long long>;
    
    int T;
    cin >> T;
    while (T--) {
        gint A, B;
        cin >> A.x >> A.y >> B.x >> B.y;
        gint res = gcd(A, B);
        cout << res.x << " " << res.y << endl;
    }
}

// New Year Contest L - 机のしみ
void New_Year_Contest_L() {
    using gint = GaussNumber<long long>;
    const long long INF = 1LL<<60;
    
    int N;
    cin >> N;
    vector<gint> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i].x >> A[i].y;
    gint g = 0;
    for (int i = 1; i < N; ++i) g = gcd(g, A[i] - A[0]);

    long long xmin = INF, ymin = INF, xmax = -INF, ymax = -INF;
    for (int i = 0; i < N; ++i) {
        gint q = A[i] / g;
        xmin = min(xmin, q.x), ymin = min(ymin, q.y);
        xmax = max(xmax, q.x), ymax = max(ymax, q.y);
    }
    long long D = max(xmax - xmin, ymax - ymin) + 1;
    cout << D * D - N << endl;
}


int main() {
    Yosupo_Gcd_of_Gaussian_Integers();
    //New_Year_Contest_L();
}

