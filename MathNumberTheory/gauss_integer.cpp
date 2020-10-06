//
// Modular Arithmetics
//
// cf.
//   https://qiita.com/drken/items/336ef288b451e86c15cb
//   https://maspypy.com/atcoder-l-%e6%9c%ba%e3%81%ae%e3%81%97%e3%81%bf%ef%bc%88new-year-contest-2015%ef%bc%89
//
// verified:
//   New Year Contest L - 机のしみ
//     https://atcoder.jp/contests/NYC2015/tasks/nyc2015_12
// 


#include <bits/stdc++.h>
using namespace std;
template<class T> inline bool chmax(T& a, T b) { if (a < b) { a = b; return 1; } return 0; }
template<class T> inline bool chmin(T& a, T b) { if (a > b) { a = b; return 1; } return 0; }


// Gauss Integer
struct gint {
    long long x, y;
    constexpr gint() : x(0), y(0) { }
    constexpr gint(long long x) : x(x), y(0) { }
    constexpr gint(long long x, long long y) : x(x), y(y) { }
    friend constexpr long long abs(const gint &r) noexcept {
        return r.x * r.x + r.y * r.y;
    }
    constexpr gint operator - () const noexcept {
        return gint(-x, -y);
    }
    constexpr gint operator + (const gint& r) const noexcept { return gint(*this) += r; }
    constexpr gint operator - (const gint& r) const noexcept { return gint(*this) -= r; }
    constexpr gint operator * (const gint& r) const noexcept { return gint(*this) *= r; }
    constexpr gint operator / (const gint& r) const noexcept { return gint(*this) /= r; }
    constexpr gint operator % (const gint& r) const noexcept { return gint(*this) %= r; }
    constexpr gint& operator += (const gint& r) noexcept {
        x += r.x, y += r.y;
        return *this;
    }
    constexpr gint& operator -= (const gint& r) noexcept {
        x -= r.x, y -= r.y;
        return *this;
    }
    constexpr gint& operator *= (const gint& r) noexcept {
        auto tx = x * r.x - y * r.y;
        auto ty = x * r.y + y * r.x;
        x = tx, y = ty;
        return *this;
    }
    constexpr gint& operator /= (const gint& r) noexcept {
        long long a = x, b = y, c = r.x, d = r.y;        
    	assert(c != 0 || d != 0);
        long long tx = (a * c + b * d) / (c * c + d * d);
        long long ty = (b * c - a * d) / (c * c + d * d);
        gint res;
        long long dmin = -1;
        for (long long nx = tx - 1; nx <= tx + 1; ++nx) {
            for (long long ny = ty - 1; ny <= ty + 1; ++ny) {
                gint q(nx, ny);
                long long d = abs((*this) - q * r);
                if (dmin == -1 || dmin > d) {
                    dmin = d;
                    res = q;
                }
            }
        }
        return *this = res;
    }
    constexpr gint& operator %= (const gint& r) noexcept {
        gint q = (*this) / r;
        return (*this) -= q * r;
    }
    constexpr bool operator == (const gint& r) const noexcept {
        return x == r.x && y == r.y;
    }
    constexpr bool operator != (const gint& r) const noexcept {
        return x != r.x || y != r.y;
    }
    friend ostream& operator << (ostream &os, const gint& r) noexcept {
        if (r.x == 0 && r.y == 0) return os << "0";
        else if (r.x == 0) return os << r.y << "i";
        else os << r.x;
        if (r.y == 0) return os;
        else if (r.y > 0) return os << " + " << r.y << "i";
        else return os << " - " << (-r.y) << "i";
    }
    friend gint gcd(const gint&x, const gint&y) {
        if (y == 0) return x;
        else return gcd(y, x % y);
    }
};


const long long INF = 1LL<<60;
int main() {
    int N;
    cin >> N;
    vector<gint> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i].x >> A[i].y;
    gint g = 0;
    for (int i = 1; i < N; ++i) g = gcd(g, A[i] - A[0]);

    long long xmin = INF, ymin = INF, xmax = -INF, ymax = -INF;
    for (int i = 0; i < N; ++i) {
        gint q = A[i] / g;
        chmin(xmin, q.x), chmin(ymin, q.y);
        chmax(xmax, q.x), chmax(ymax, q.y);
    }
    long long D = max(xmax - xmin, ymax - ymin) + 1;
    cout << D * D - N << endl;
}
