//
// floor sum
//
// verified
//   AtCoder ABC 283 Ex - Popcount Sum
//     https://atcoder.jp/contests/abc283/tasks/abc283_h
//
//   AOJ 3217 Cafe au lait (OUPC 2020 I)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=3217
//
//   yukicoder No.2066 Simple Math !
//     https://yukicoder.me/problems/no/2066
//


#include <bits/stdc++.h>
using namespace std;


// sum_{i=0}^{n-1} floor((a * i + b) / m)
// O(log(n + m + a + b))
// __int128 can be used for T
template<class T> T floor_sum(T n, T a, T b, T m) {
    if (n == 0) return 0;
    T res = 0;
    if (a >= m) {
        res += n * (n - 1) * (a / m) / 2;
        a %= m;
    }
    if (b >= m) {
        res += n * (b / m);
        b %= m;
    }
    if (a == 0) return res;
    T ymax = (a * n + b) / m, xmax = ymax * m - b;
    if (ymax == 0) return res;
    res += (n - (xmax + a - 1) / a) * ymax;
    res += floor_sum(ymax, m, (a - xmax % a) % a, a);
    return res;
}

// #lp under (and on) the segment (x1, y1)-(x2, y2)
// not including y = 0, x = x2
template<class T> T num_lattice_points(T x1, T y1, T x2, T y2) {
    T dx = x2 - x1;
    return floor_sum(dx, y2 - y1, dx * y1, dx);
}


/////////////////////////////////////////
// Solvers
/////////////////////////////////////////

/* yukicoder No.2066 */

// calc #n that can be expressed n =  Px + Qy (P, Q is coprime)
// 0 <= n <= M
long long calc_num(__int128 P, __int128 Q, __int128 M) {
    __int128 mp = M / P;
    __int128 N = min(mp + 1, Q);
    __int128 a = P, b = M + Q - a * (N - 1);
    return floor_sum(N, a, b, Q) - 1;
}

void solveYukicoder2066() {
    int CASE;
    cin >> CASE;
    while (CASE--) {
        long long P, Q, K;
        cin >> P >> Q >> K;

        long long G = gcd(P, Q);
        P /= G, Q /= G;

        long long low = -1, high = 1LL<<50;
        while (high - low > 1) {
            long long M = (low + high) / 2;
            if (calc_num(P, Q, M) >= K) high = M;
            else low = M;
        }
        cout << high * G << endl;
    }
}


/* AOJ 3217 */
// modint
template<int MOD> struct Fp {
    long long val;
    constexpr Fp(long long v = 0) noexcept : val(v % MOD) {
        if (val < 0) val += MOD;
    }
    constexpr int getmod() const { return MOD; }
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
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
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
    friend constexpr istream& operator >> (istream& is, Fp<MOD>& x) noexcept {
        is >> x.val;
        x.val %= MOD;
        if (x.val < 0) x.val += MOD;
        return is;
    }
    friend constexpr ostream& operator << (ostream& os, const Fp<MOD>& x) noexcept {
        return os << x.val;
    }
    friend constexpr Fp<MOD> modpow(const Fp<MOD>& r, long long n) noexcept {
        if (n == 0) return 1;
        if (n < 0) return modpow(modinv(r), -n);
        auto t = modpow(r, n / 2);
        t = t * t;
        if (n & 1) t = t * r;
        return t;
    }
    friend constexpr Fp<MOD> modinv(const Fp<MOD>& r) noexcept {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        return Fp<MOD>(u);
    }
};

void solveAOJ3217() {
    const int MOD = 1000000007;
    using mint = Fp<MOD>;

    int N;
    cin >> N;
    vector<long long> X(N), Y(N);
    for (int i = 0; i < N; ++i) cin >> X[i] >> Y[i];
    vector<int> ids(N);
    iota(ids.begin(), ids.end(), 0);
    sort(ids.begin(), ids.end(), [&](int i, int j) {
        return Y[i]*X[j] < Y[j]*X[i];
    });
    mint res = 0;
    long long sy = 0;
    for (auto i : ids) {
        res -= mint(X[i]) * mint(sy);
        res -= num_lattice_points(0LL, 0LL, X[i], Y[i]);
        res += gcd(X[i], Y[i]);
        sy += Y[i];
    }
    reverse(ids.begin(), ids.end());
    sy = 0;
    for (auto i : ids) {
        res += mint(X[i]) * mint(sy);
        res += num_lattice_points(0LL, 0LL, X[i], Y[i]);
        sy += Y[i];
    }
    cout << res << endl;
}

int main() {
    //solveYukicoder2066();
    solveAOJ3217();
}