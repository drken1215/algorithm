#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;

// chmax, chmin
template<class T> inline bool chmax(T& a, T b) { if (a < b) { a = b; return 1; } return 0; }
template<class T> inline bool chmin(T& a, T b) { if (a > b) { a = b; return 1; } return 0; }

// debug stream of pair, vector 
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }


// 有理数
long long calc_gcd(long long a, long long b) {return b ? calc_gcd(b, a % b) : a;}
struct frac {
    long long first, second;

    using D = long double;
    inline frac normalize() {
        if (second < 0) {first = -first; second = -second;}
        long long d = calc_gcd(abs(first), abs(second));
        if (d == 0) {first = 0; second = 1;}
        else {first /= d; second /= d;}
        return *this;
    }
    frac(long long f = 0, long long s = 1) : first(f), second(s) { normalize(); }
    inline D to_d() const { return D(first) / second; }
    inline frac operator - () { (*this).first *= -1; return (*this); }
    inline const frac& operator = (long long a) { *this = frac(a, 1); return *this; }
    inline const frac& operator += (const frac& a);
    inline const frac& operator += (long long a);
    inline const frac& operator -= (const frac& a);
    inline const frac& operator -= (long long a);
    inline const frac& operator *= (const frac& a);
    inline const frac& operator *= (long long a);
    inline const frac& operator /= (const frac& a);
    inline const frac& operator /= (long long a);
    inline friend ostream& operator << (ostream& s, const frac& f) { 
        s << f.first; if (f.second != 1) s << "/" << f.second; return s;
    }
};
inline bool operator == (const frac &a, const frac&b) {
    return a.first * b.second == a.second * b.first;
}
inline bool operator != (const frac &a, const frac &b) { return !(a == b); }
inline bool operator < (const frac& a, const frac& b) {
    return a.first * b.second < a.second * b.first;
}
inline bool operator > (const frac& a, const frac& b) { return b < a; }
inline bool operator <= (const frac& a, const frac& b) {
    return a.first * b.second <= a.second * b.first;
}
inline bool operator >= (const frac& a, const frac& b) { return b <= a; }
inline frac operator + (const frac& a, const frac& b) {
    frac res;
    res.first = a.first * b.second + a.second * b.first;
    res.second = a.second * b.second;
    res.normalize();
    return res;
}
inline frac operator - (const frac& a, const frac& b) {
    frac res;
    res.first = a.first * b.second - a.second * b.first;
    res.second = a.second * b.second;
    res.normalize();
    return res;
}
inline frac operator * (const frac& a, const frac& b) {
    frac res;
    res.first = a.first * b.first;
    res.second = a.second * b.second;
    res.normalize();
    return res;
}
inline frac operator / (const frac& a, const frac& b) {
    frac res;
    res.first = a.first * b.second;
    res.second = a.second * b.first;
    res.normalize();
    return res;
}
inline frac abs(const frac& a) {
    frac res; res = a; res.normalize(); 
    if (res.first < 0) res.first = res.first * (-1);
    return res;
}
inline const frac& frac::operator += (const frac& x) {*this = *this + x; return *this;}
inline const frac& frac::operator += (long long x) {*this = *this + x; return *this;}
inline const frac& frac::operator -= (const frac& x) {*this = *this - x; return *this;}
inline const frac& frac::operator -= (long long x) {*this = *this + x; return *this;}
inline const frac& frac::operator *= (const frac& x) {*this = *this * x; return *this;}
inline const frac& frac::operator *= (long long x) {*this = *this * x; return *this;}
inline const frac& frac::operator /= (const frac& x) {*this = *this / x; return *this;}
inline const frac& frac::operator /= (long long x) {*this = *this / x; return *this;}



long double cost(const frac &f, long long T) {
    long double df = f.to_d();
    return sqrt(df * T * T / 2 + 0.5 / df);
}

int main() {
    long long T; int N;
    cin >> T >> N;
    const frac center(1, T); // 45 度打ち出しの場合
    using pf = pair<frac,frac>;
    vector<pf> upper, lower, middle;
    for (int i = 0; i < N; ++i) {
        long long x, low, up; cin >> x >> low >> up;
        frac flow(low, x * (T-x));
        frac fup(up, x * (T-x));
        if (fup < center) lower.push_back({flow, fup});
        else if (flow > center) upper.push_back({flow, fup});
        else middle.push_back({flow, fup});
    }
    sort(lower.begin(), lower.end(), [](const pf &a, const pf &b) {
            return a.second < b.second;});
    sort(upper.begin(), upper.end(), [](const pf &a, const pf &b) {
            return a.first > b.first;});

    long double res = 0.0;
    frac left = 0, right = 1000100; // right * 10^12/4 がオーバーフローしないように
    for (auto inter : lower) {
        if (left >= inter.first) continue;
        res += cost(inter.second, T);
        left = inter.second;
    }
    for (auto inter : upper) {
        if (right <= inter.second) continue;
        res += cost(inter.first, T);
        right = inter.first;
    }
    bool remain = false;
    for (auto inter : middle) {
        if (left < inter.first && inter.second < right) remain = true;
    }
    if (remain) res += cost(center, T);
    cout << fixed << setprecision(20) << res << endl;
}
