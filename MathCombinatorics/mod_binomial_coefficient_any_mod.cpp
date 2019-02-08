//
// nCr mod m (m: 任意, n, r <= 10^6)
//
// verified:
//   ARC 012 D - Don't worry. Be Together
//     https://atcoder.jp/contests/arc012/tasks/arc012_4
//


#include <iostream>
#include <vector>
#include <map>
#include <cstring>
using namespace std;


const int MAX = 110000;

// p is prime, m = p^k
// a! = p^(ord[a]) * fac[a] (mod. m)
long long ord[MAX], fac[MAX];
void prime_com_init(long long p, long long pm) {
    ord[0] = ord[1] = 0;
    fac[0] = fac[1] = 1;
    for (int i = 2; i < MAX; i++) {
        long long add = 0;
        long long ni = i;
        while (ni % p == 0) ++add, ni /= p;
        ord[i] = ord[i-1] + add;
        fac[i] = fac[ni-1] * ni % pm;
    }
}

inline long long mod(long long a, long long m) {
    return (a % m + m) % m;
}

long long modpow(long long a, long long n, long long mod) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}

long long modinv(long long a, long long mod) {
    long long b = mod, u = 1, v = 0;
    while (b) {
        long long t = a/b;
        a -= t*b; swap(a, b);
        u -= t*v; swap(u, v);
    }
    u %= mod;
    if (u < 0) u += mod;
    return u;
}

// nCr mod. pm
long long COM(long long n, long long r, long long p, long long pm) {
    if (n < 0 || r < 0 || n < r) return 0;
    long long e = ord[n] - ord[r] - ord[n-r];
    long long res = fac[n] * modinv(fac[r] * fac[n-r] % pm, pm) % pm;
    res = res * modpow(p, e, pm) % pm;
    return res;
}

long long extGcd(long long a, long long b, long long &p, long long &q) {
    if (b == 0) { p = 1; q = 0; return a; }
    long long d = extGcd(b, a%b, q, p);
    q -= a/b * p;
    return d;
}

long long Garner(vector<long long> b, vector<long long> m, long long MOD) {
    m.push_back(MOD); // banpei
    vector<long long> coeffs((int)m.size(), 1);
    vector<long long> constants((int)m.size(), 0);
    for (int k = 0; k < (int)b.size(); ++k) {
        long long t = mod((b[k] - constants[k]) * modinv(coeffs[k], m[k]), m[k]);
        for (int i = k+1; i < (int)m.size(); ++i) {
            (constants[i] += t * coeffs[i]) %= m[i];
            (coeffs[i] *= m[k]) %= m[i];
        }
    }
    return constants.back();
}

vector<pair<long long, long long> > prime_factorize(long long n) {
    vector<pair<long long, long long> > res;
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p != 0) continue;
        int num = 0;
        while (n % p == 0) { ++num; n /= p; }
        res.push_back(make_pair(p, num));
    }
    if (n != 1) res.push_back(make_pair(n, 1));
    return res;
}



int N, T, M;
vector<int> x, y, r;

long long solve() {
    cin >> N >> T >> M;
    x.resize(N); y.resize(N); r.resize(N);
    bool ok = true;
    for (int i = 0; i < N; ++i) {
        cin >> x[i] >> y[i];
        if (x[i] < 0) x[i] = -x[i];
        if (y[i] < 0) y[i] = -y[i];
        if (T < x[i] + y[i]) ok = false;
        if ((T - x[i] - y[i]) % 2 != 0) ok = false;
        r[i] = (T - x[i] - y[i]) / 2;
    }
    if (!ok) return 0;
    if (M == 1) return 0;
    
    vector<pair<long long, long long> > pf = prime_factorize(M);
    vector<long long> vb, vm;
    for (auto ps : pf) {
        long long p = ps.first, e = ps.second;
        long long pm = 1;
        for (int i = 0; i < e; ++i) pm *= p;
        prime_com_init(p, pm);
        long long b = 1;
        for (int i = 0; i < N; ++i) {
            b *= COM(T, r[i], p, pm) * COM(T, x[i] + r[i], p, pm) % pm;
            b %= pm;
        }
        vm.push_back(pm);
        vb.push_back(b);
    }
    auto res = Garner(vb, vm, M);
    return res;
}

int main() {
    cout << solve() << endl;
}
