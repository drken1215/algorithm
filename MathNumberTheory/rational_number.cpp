//
// 有理数
//
// verified
//   AOJ 1426 Remodeling the Dungeon 2 (ICPC Asia 2024 H)
//     https://onlinejudge.u-aizu.ac.jp/problems/1462
//


#include <bits/stdc++.h>
using namespace std;


// 有理数
template<class T> struct frac {
    // inner values
    T first, second;

    // constructor
    frac& normalize() {
        if (second < 0) first = -first, second = -second;
        T d = gcd(abs(first), abs(second));
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


int main() {
    ICPC_Asia_2024_J();
}


// int main() {
//     long long T; int N;
//     cin >> T >> N;
//     const frac center(1, T); // 45 度打ち出しの場合
//     using pf = pair<frac,frac>;
//     vector<pf> upper, lower, middle;
//     for (int i = 0; i < N; ++i) {
//         long long x, low, up; cin >> x >> low >> up;
//         frac flow(low, x * (T-x));
//         frac fup(up, x * (T-x));
//         if (fup < center) lower.push_back({flow, fup});
//         else if (flow > center) upper.push_back({flow, fup});
//         else middle.push_back({flow, fup});
//     }
//     sort(lower.begin(), lower.end(), [](const pf &a, const pf &b) {
//             return a.second < b.second;});
//     sort(upper.begin(), upper.end(), [](const pf &a, const pf &b) {
//             return a.first > b.first;});

//     long double res = 0.0;
//     frac left = 0, right = 1000100; // right * 10^12/4 がオーバーフローしないように
//     for (auto inter : lower) {
//         if (left >= inter.first) continue;
//         res += cost(inter.second, T);
//         left = inter.second;
//     }
//     for (auto inter : upper) {
//         if (right <= inter.second) continue;
//         res += cost(inter.first, T);
//         right = inter.first;
//     }
//     bool remain = false;
//     for (auto inter : middle) {
//         if (left < inter.first && inter.second < right) remain = true;
//     }
//     if (remain) res += cost(center, T);
//     cout << fixed << setprecision(20) << res << endl;
// }