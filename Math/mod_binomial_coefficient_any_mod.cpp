//
// nCr mod p^m (p: 素数, p^m <= 10^9,  n, r <= 10^6)
//   中国剰余定理と組み合わせることで、任意 mod nCr も可能 
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


// mod function
long long mod(long long a, long long mod) {
    return (a % mod + mod) % mod;
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
        a -= t * b, swap(a, b);
        u -= t * v, swap(u, v);
    }
    u %= mod;
    if (u < 0) u += mod;
    return u;
}

// binomial coefficient
struct BiCoef {
    // max size of n of nCr
    int n_;
    
    // pm = p^k, p is prime
    // mod pm
    long long p_, pm_;
    
    // i! = p^ord[i] * fact[i] (mod. m)
    vector<long long> ord_, fact_;

    // constructor
    BiCoef(int n) : n_(n), ord_(n), fact_(n) {}
    BiCoef(long long p, long long pm, int n) :
        n_(n), p_(p), pm_(pm), ord_(n), fact_(n) {
        init(p, pm);
    }
    void init(int n) {
        ord_.resize(n);
        fact_.resize(n);
    }
    void init(long long p, long long pm) {
        p_ = p, pm_ = pm;
        ord_[0] = ord_[1] = 0;
        fact_[0] = fact_[1] = 1;
        for (int i = 2; i < n_; i++) {
            long long add = 0;
            long long ni = i;
            while (ni % p == 0) ++add, ni /= p;
            ord_[i] = ord_[i-1] + add;
            fact_[i] = fact_[ni-1] * ni % pm;
        }
    }
    void init(long long p, long long pm, int n) {
        init(n);
        init(p, pm);
    }

    // nCr mod. pm
    long long com(long long n, long long r) {
        if (n < 0 || r < 0 || n < r) return 0;
        long long e = ord_[n] - ord_[r] - ord_[n-r];
        long long res = fact_[n] * modinv(fact_[r] * fact_[n-r] % pm_, pm_) % pm_;
        res = res * modpow(p_, e, pm_) % pm_;
        return res;
    }
};

// Garner's Algorithm
long long Garner(vector<long long> b, vector<long long> m, long long MOD) {
    m.push_back(MOD); // banpei
    vector<long long> coeffs(m.size(), 1);
    vector<long long> constants(m.size(), 0);
    for (int k = 0; k < (int)b.size(); ++k) {
        long long t = mod((b[k] - constants[k]) * modinv(coeffs[k], m[k]), m[k]);
        for (int i = k+1; i < (int)m.size(); ++i) {
            (constants[i] += t * coeffs[i]) %= m[i];
            (coeffs[i] *= m[k]) %= m[i];
        }
    }
    return constants.back();
}



//------------------------------//
// Examples
//------------------------------//

// 入力
int N, T, M;
vector<int> x, y, r;

// solver
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

    // 素因数分解
    vector<long long> ps, pms;
    long long oriM = M;
    for (long long p = 2; p * p <= M; ++p) {
        if (M % p != 0) continue;
        long long pm = 1;
        while (M % p == 0) {
            pm *= p;
            M /= p;
        }
        ps.push_back(p), pms.push_back(pm);
    }
    if (M != 1) ps.push_back(M), pms.push_back(M);
    
    // nCr mod pm
    BiCoef bf(110000);
    vector<long long> bs;
    for (int i = 0; i < ps.size(); ++i) {
        long long p = ps[i], pm = pms[i];
        bf.init(p, pm);
        long long b = 1;
        for (int i = 0; i < N; ++i) {
            b *= bf.com(T, r[i]) * bf.com(T, x[i] + r[i]) % pm;
            b %= pm;
        }
        bs.push_back(b);
    }
    auto res = Garner(bs, pms, oriM);
    return res;
}

int main() {
    cout << solve() << endl;
}
