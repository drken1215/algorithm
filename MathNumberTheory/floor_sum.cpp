//
// floor sum
//
// verified
//   OUPC 2020 I - カフェオレ
//     https://onlinejudge.u-aizu.ac.jp/beta/room.html#OUPC2020/problems/I
//


#include <bits/stdc++.h>
using namespace std;


// sum_{i=0}^{n-1} floor((a * i + b) / m)
// O(log(n + m + a + b))
long long floor_sum(long long n, long long a, long long b, long long m) {
    if (n == 0) return 0;
    long long res = 0;
    if (a >= m) {
        res += n * (n - 1) * (a / m) / 2;
        a %= m;
    }
    if (b >= m) {
        res += n * (b / m);
        b %= m;
    }
    if (a == 0) return res;
    long long ymax = (a * n + b) / m, xmax = ymax * m - b;
    if (ymax == 0) return res;
    res += (n - (xmax + a - 1) / a) * ymax;
    res += floor_sum(ymax, m, (a - xmax % a) % a, a);
    return res;
}

// #lp under (and on) the segment (x1, y1)-(x2, y2)
// not including y = 0, x = x2
long long num_lattice_points(long long x1, long long y1, long long x2, long long y2) {
    long long dx = x2 - x1;
    return floor_sum(dx, y2 - y1, dx * y1, dx);
}


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

const int MOD = 1000000007;
using mint = Fp<MOD>;

long long GCD(long long x, long long y) {
    return y ? GCD(y, x % y) : x;
}

int main() {
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
        res -= num_lattice_points(0, 0, X[i], Y[i]);
        res += GCD(X[i], Y[i]);
        sy += Y[i];
    }
    reverse(ids.begin(), ids.end());
    sy = 0;
    for (auto i : ids) {
        res += mint(X[i]) * mint(sy);
        res += num_lattice_points(0, 0, X[i], Y[i]);
        sy += Y[i];
    }
    cout << res << endl;
}
