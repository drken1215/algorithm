//
// Pollard のロー素因数分解法
//
// cf:
//   素因数分解を O(n^(1/4)) でする (Kiri)
//     https://qiita.com/Kiri8128/items/eca965fe86ea5f4cbb98
//
// verifed:
//   Yosupo Judge Factorize
//     https://judge.yosupo.jp/problem/factorize
//
//   アルゴ式 番外編：ポラードのロー素因数分解法
//     https://algo-method.com/tasks/553
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
// xor128 rng
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


//-///////////////////////////-//
// Example
//-///////////////////////////-//

void YosupoJudgeFactorize() {
    // 入力
    int N;
    cin >> N;
    // 素因数分解
    for (int i = 0; i < N; ++i) {
        long long a;
        cin >> a;
        const auto& res = prime_factorize(a);
        cout << res.size();
        for (auto p : res) cout << " " << p;
        cout << endl;
    }
}


int main() {
    YosupoJudgeFactorize();
}

