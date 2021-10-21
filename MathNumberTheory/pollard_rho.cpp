//
// Pollard のロー素因数分解法
//
// verifed
//   アルゴ式 番外編：ポラードのロー素因数分解法
//     https://algo-method.com/tasks/553
//     


#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;


// Miller-Rabin 素数判定法
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

// Pollard のロー法
long long gcd(long long A, long long B) {
    A = abs(A), B = abs(B);
    if (B == 0) return A;
    else return gcd(B, A % B);
}
    
long long pollard(long long N) {
    if (N % 2 == 0) return 2;
    if (is_prime(N)) return N;

    auto f = [&](long long x) -> long long {
        return (__int128_t(x) * x + 1) % N;
    };
    long long step = 0;
    while (true) {
        ++step;
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


int main() {
    // 入力
    int N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];

    // 素因数分解
    for (auto a : A) {
        const auto& res = prime_factorize(a);
        for (auto p : res) cout << p << " ";
        cout << endl;
    }
}
