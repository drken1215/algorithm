//
// 確率的素数判定 (Miller-Rabin)
//
// verifed
//   [アルゴ式 番外編：ミラー-ラビン素数判定法](https://algo-method.com/tasks/513)
//     


#include <iostream>
using namespace std;


// A * B を M で割ったあまり
long long mod_mul(long long A, long long B, long long M) {
    if (B == 0) return 0;
    long long res = mod_mul((A + A) % M, B >> 1, M);
    if (B & 1) res = (res + A) % M;
    return res;
}

// A ^ N を M で割ったあまり
long long mod_pow(long long A, long long N, long long M) {
    if (N == 0) return 1 % M;
    long long res = mod_pow(A, N / 2, M);
    res = mod_mul(res, res, M);
    if (N & 1) res = mod_mul(res, A, M);
    return res;
}

// ミラー・ラビン素数判定法
// N < 341,550,071,728,321 (17 まで) については確定的に判定可能
// N < 3,825,123,056,546,413,051 (23 まで) についても確定的に判定可能
bool is_prime(long long N) {
    if (N == 2) return true;
    if (!(N & 1) || N <= 1) return false;
    vector<long long> A = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    long long s = 0, d = N - 1;
    while (!(d & 1)) {
        ++s;
        d = d >> 1;
    }
    for (int i = 0; i < A.size() && A[i] < N; ++i) {
        long long x = mod_pow(A[i], d, N), t;
        if (x != 1) {
            for (t = 0; t < s; ++t) {
                if (x == N - 1) break;
                x = mod_mul(x, x, N);
            }
            if (t == s) return false;
        }
    }
    return true;
}

int main() {
    // 入力
    int N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];

    // 素数判定
    for (auto a : A) {
        if (is_prime(a))
            cout << "Yes" << endl;
        else
            cout << "No" << endl;
    }
}
