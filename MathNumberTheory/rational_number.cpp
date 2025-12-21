//
// 有理数
//
// verified
//   AOJ 1426 Remodeling the Dungeon 2 (ICPC Asia 2024 H)
//     https://onlinejudge.u-aizu.ac.jp/problems/1464
//
//   数学アルゴ本 071 - Linear Programming
//     https://atcoder.jp/contests/math-and-algorithm/tasks/math_and_algorithm_bf
//


#include <bits/stdc++.h>
using namespace std;


// 有理数
template<class T> struct frac {
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
        T abs_first = (first >= 0 ? first : -first);
        T d = gcd(abs_first, second);
        if (d == 0) first = 0, second = 1;
        else first /= d, second /= d;
        return *this;
    }
    constexpr frac(T f = 0, T s = 1) : first(f), second(s) { 
        normalize(); 
    }
    constexpr frac& operator = (T a) { 
        *this = frac(a, 1); 
        return *this;
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

// AOJ 1426 Remodeling the Dungeon 2 (ICPC Asia 2024 H)
void ICPC_Asia_2024_J() {
    using lfrac = frac<long long>;
    const long long M = 10000;
    long long N, S, C, sumA = 0;
    cin >> N >> S >> C;
    vector<long long> A(N), L(N), R(N);
    for (int i = 0; i < N; i++) cin >> A[i] >> L[i] >> R[i], sumA += A[i];

    vector<lfrac> vw({lfrac(-1), lfrac(1)});
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            lfrac f(-L[i]+R[i]+L[j]-R[j], L[i]+R[i]-L[j]-R[j]);
            if (f >= -1 && f <= 1) vw.push_back(f);
        }
    }
    sort(vw.begin(), vw.end());
    vw.erase(unique(vw.begin(), vw.end()), vw.end());

    auto func = [&](lfrac w) -> lfrac { 
        vector<pair<lfrac,int>> v;
        for (int i = 0; i < N; i++) {
            v.emplace_back(w*lfrac(L[i]+R[i])/2 + lfrac(L[i]-R[i])/2, i);
        }
        sort(v.begin(), v.end());

        vector<lfrac> y(N, -1);
        lfrac z, limit = sumA - S;
        for (auto [val, i] : v) {
            if (limit >= 0) {
                limit -= A[i];
                y[i] = 0;
                if (limit < 0) z = -val;
            } else {
                y[i] = z + val;
            }
        }
        lfrac obj = z * S + w * C * S;
        for (int i = 0; i < N; i++) obj -= y[i] * A[i];
        return obj;
    };

    int left = 0, right = vw.size()-1;
    while (right - left > 2) {
        int m1 = (left * 2 + right) / 3, m2 = (left + right * 2) / 3;
        if (func(vw[m1]) > func(vw[m2])) right = m2;
        else left = m1;
    }
    lfrac res = max({func(vw[left]), func(vw[(left+right)/2]), func(vw[right])});
    res /= M;
    cout << res.first << " " << res.second << endl;
}

// 数学アルゴ本 071 - Linear Programming
using i128 = __int128;
constexpr i128 to_integer(const string &s) {
    i128 res = 0;
    for (auto c : s) {
         if (isdigit(c)) res = res * 10 + (c - '0');
    }
    if (s[0] == '-') res *= -1;
    return res;
}
istream& operator >> (istream &is, i128 &x) {
    string s;
    is >> s;
    x = to_integer(s);
    return is;
}
void MathAlgorithm071() {
    using FF = frac<i128>;
    int N;
    cin >> N;
    vector<i128> A(N), B(N), C(N);
    for (int i = 0; i < N; i++) cin >> A[i] >> B[i] >> C[i];

    auto check = [&](FF x, FF y) {
        for (int i = 0; i < N; i++) {
            if (x * A[i] + y * B[i] > C[i]) return false;
        }
        return true;
    };

    FF res = i128(0);
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            FF a = A[i], b = B[i], c = A[j], d = B[j], k = C[i], l = C[j];

            if (a * d - b * c == i128(0)) continue;
            FF det = a * d - b * c;
            FF x = (d * k - b * l) / det, y = (a * l - c * k) / det;
            if (check(x, y)) res = max(res, x + y);
        }
    }
    using DD = long double;
    DD dres = DD(res.first) / res.second;
    cout << fixed << setprecision(10) << dres << endl;
}


int main() {
    //ICPC_Asia_2024_J();
    MathAlgorithm071();
}