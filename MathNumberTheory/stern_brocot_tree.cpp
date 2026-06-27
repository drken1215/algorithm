//
// Stern-Brocot 木上の二分探索
//
// verified
//   Library Checker - Rational Approximation
//     https://judge.yosupo.jp/problem/rational_approximation
//
//   ICPC アジア地区 京都大会 1999 A - Rational Irrationals (AOJ 1208)
//     https://onlinejudge.u-aizu.ac.jp/problems/1208
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


// Stern-Brocot Tree
template<class T> struct SternBrocotTree {
    // binary search on Stern-Brocot Tree
    // return {l (= a/b), r (= c/d(} s.t. l: OK, r: NG
    // and a, b, c, d are maximized where a, b, c, d <= lim
    template<class Func> static tuple<T, T, T, T> binary_search(Func check, T lim) {
        assert(check(0, 1));
        assert(!check(1, 0));
        auto rec = [&](auto &&rec, bool which, T &a, T &b, T c, T d) -> void {
            if (a + c > lim || b + d > lim) return;
            if (check(a + c, b + d) == which) {
                a += c, b += d;
                rec(rec, which, a, b, c + c, d + d);
            }
            if (a + c <= lim && b + d <= lim && check(a + c, b + d) == which) a += c, b += d;
        };
        T a = 0, b = 1, c = 1, d = 0;
        while (a + c <= lim && b + d <= lim) {
            rec(rec, true, a, b, c, d);
            rec(rec, false, c, d, a, b);
        }
        return {a, b, c, d};
    }
};

// not necessary, but often use: rational number
template<class T = long long> struct frac {
    // gcd
    static T gcd(T a, T b) {
        a = max(a, -a), b = max(b, -b);
        while (b) {
            a %= b;
            swap(a, b);
        }
        return a;
    }

    // inner values
    T first, second;

    // constructor
    frac& normalize() {
        if (first == 0 && second != 0) {
            second = 1;
            return *this;
        }
        if (second == 0 && first != 0) {
            first = 1;
            return *this;
        }
        if (second < 0) first = -first, second = -second;
        T d = gcd(max(first, -first), second);
        if (d == 0) first = 0, second = 1;
        else first /= d, second /= d;
        return *this;
    }
    frac(const frac&) = default;
    frac& operator = (const frac&) = default;
    constexpr frac(T f = 0, T s = 1) : first(f), second(s) { 
        normalize(); 
    }
    constexpr frac& operator = (T a) { 
        *this = frac(a, 1); 
        return *this;
    }
    constexpr long double to_double() const {
        assert(second != 0);
        return (long double)(first) / (long double)(second);
    }
    friend constexpr long double to_double(const frac &r) {
        return r.to_double();
    }

    // comparison operators
    constexpr bool operator == (const frac &r) const {
        return this->first == r.first && this->second == r.second;
    }
    constexpr bool operator != (const frac &r) const {
        return this->first != r.first || this->second != r.second;
    }
    constexpr bool operator < (const frac &r) const {
        return this->first * r.second < this->second * r.first;
    }
    constexpr bool operator > (const frac &r) const {
        return this->first * r.second > this->second * r.first;
    }
    constexpr bool operator <= (const frac &r) const {
        return this->first * r.second <= this->second * r.first;
    }
    constexpr bool operator >= (const frac &r) const {
        return this->first * r.second >= this->second * r.first;
    }
    
    // arithmetic operators
    constexpr frac& operator += (const frac &r) {
        this->first = this->first * r.second + this->second * r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator -= (const frac &r) {
        this->first = this->first * r.second - this->second * r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator *= (const frac &r) {
        this->first *= r.first;
        this->second *= r.second;
        this->normalize();
        return *this;
    }
    constexpr frac& operator /= (const frac &r) {
        this->first *= r.second;
        this->second *= r.first;
        this->normalize();
        return *this;
    }
    constexpr frac operator + () const { return frac(*this); }
    constexpr frac operator - () const { return frac(0) - frac(*this); }
    constexpr frac operator + (const frac &r) const { return frac(*this) += r; }
    constexpr frac operator - (const frac &r) const { return frac(*this) -= r; }
    constexpr frac operator * (const frac &r) const { return frac(*this) *= r; }
    constexpr frac operator / (const frac &r) const { return frac(*this) /= r; }
    friend constexpr ostream& operator << (ostream &os, const frac<T> &x) {
        os << x.first; 
        if (x.second != 1) os << "/" << x.second;
        return os;
    }
};


//------------------------------//
// Examples
//------------------------------//

// Library Checker - Rational Approximation
void LibraryCheckerRatilnalApproximation() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    using sbt = SternBrocotTree<long long>;
    int T;
    cin >> T;
    while (T--) {
        long long N, x, y;
        cin >> N >> x >> y;
        auto check = [&](long long a, long long b) -> bool { return a * y <= b * x; };
        auto [a, b, c, d] = sbt::binary_search(check, N);
        if (a * y == b * x) cout << a << ' ' << b << ' ' << a << ' ' << b << '\n';
        else cout << a << ' ' << b << ' ' << c << ' ' << d << '\n';
    }
}

// ICPC アジア地区 京都大会 1999 A - Rational Irrationals (AOJ 1208)
void AOJ_1208() {
    using sbt = SternBrocotTree<long long>;
    long long P, N;
    while (cin >> P >> N, P) {
        auto check = [&](long long a, long long b) -> bool {
            return a*a < b*b*P;
        };
        auto [a, b, c, d] = sbt::binary_search(check, N);
        cout << c << "/" << d << " " << a << "/" << b << endl;
    }
}


int main() {
    //LibraryCheckerRatilnalApproximation();
    AOJ_1208();
}