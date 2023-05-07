//
// (最小の) 原始根を求める
//
// verified:
//   Yosupo Judge Primitive Root
//     https://judge.yosupo.jp/problem/primitive_root
//


#include <bits/stdc++.h>
using namespace std;

// Miller-Rabin
template<class T> T pow_mod(T A, T N, T M) {
    T res = 1 % M;
    A %= M;
    while (N) {
        if (N & 1) res = (res * A) % M;
        A = (A * A) % M;
        N >>= 1;
    }
    return res;
}

bool is_prime(long long N) {
    if (N <= 1) return false;
    if (N == 2 || N == 3) return true;
    if (N % 2 == 0) return false;
    vector<long long> A = {2, 325, 9375, 28178, 450775,
                           9780504, 1795265022};
    long long s = 0, d = N - 1;
    while (d % 2 == 0) {
        ++s;
        d >>= 1;
    }
    for (auto a : A) {
        if (a % N == 0) return true;
        long long t, x = pow_mod<__int128_t>(a, d, N);
        if (x != 1) {
            for (t = 0; t < s; ++t) {
                if (x == N - 1) break;
                x = __int128_t(x) * x % N;
            }
            if (t == s) return false;
        }
    }
    return true;
}

// Pollard's Rho
unsigned int xor_shift_rng() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty; ty = tz; tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}

int xor_shift_rng(int minv, int maxv) {
    return xor_shift_rng() % (maxv - minv + 1) + minv;
}

long long pollard(long long N) {
    if (N % 2 == 0) return 2;
    if (is_prime(N)) return N;
    
    long long r = xor_shift_rng();  // random r
    auto f = [&](long long x) -> long long {
        return (__int128_t(x) * x + r) % N;
    };
    long long step = 0;
    while (true) {
        ++step;
        r = xor_shift_rng();
        long long x = step, y = f(x);
        while (true) {
            long long p = gcd(y - x + N, N);
            if (p == 0 || p == N) break;
            if (p != 1) return p;
            x = f(x);
            y = f(f(y));
        }
    }
}

vector<long long> prime_factorize(long long N) {
    if (N == 1) return {};
    long long p = pollard(N);
    if (p == N) return {p};
    vector<long long> left = prime_factorize(p);
    vector<long long> right = prime_factorize(N / p);
    left.insert(left.end(), right.begin(), right.end());
    sort(left.begin(), left.end());
    return left;
}

long long calc_primitive_root(long long p) {
    if (p == 1) return -1;
    if (p == 2) return 1;
    if (p == 998244353) return 3;
    
    const auto &pftmp = prime_factorize(p - 1);
    vector<long long> pf;
    for (auto q : pftmp) {
        if (pf.empty() || pf.back() != q) pf.push_back(q);
    }
    for (long long g = 1; g < p; g++) {
        bool ok = true;
        for (auto q : pf) {
            if (pow_mod<__int128_t>(g, (p - 1) / q, p) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
    return -1;
}


//-///////////////////////////-//
// Example
//-///////////////////////////-//

int main() {
    int Q;
    cin >> Q;
    for (int i = 0; i < Q; ++i) {
        long long P;
        cin >> P;
        cout << calc_primitive_root(P) << endl;
    }
}


