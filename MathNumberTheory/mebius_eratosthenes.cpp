//
// エラトスネテスの篩を用いて、1〜n のメビウス関数を高速に求める
//
// cf.
//   高速素因数分解
//     http://www.osak.jp/diary/diary_201310.html#20131017
//
// verified
//   Codeforces #548 Div. 2 D - Steps to One
//     https://codeforces.com/contest/1139/problem/D
//


#include <iostream>
#include <vector>
using namespace std;


// エラトステネスやりながら、「最小の素因数」「メビウス関数」を求める
const int MAX = 101010;
bool is_prime[MAX];
int min_factor[MAX];
int meb[MAX];
vector<int> preprocess(int n = MAX) {
    vector<int> res;
    for (int i = 0; i < n; ++i) is_prime[i] = true, min_factor[i] = -1,  meb[i] = 1;
    is_prime[0] = false; is_prime[1] = false;
    min_factor[0] = 0, min_factor[1] = 1;
    meb[0] = 0; meb[1] = 1;
    for (int i = 2; i < n; ++i) {
        if (!is_prime[i]) continue;
        res.push_back(i);
        min_factor[i] = i;
        meb[i] = -1;
        for (int j = i*2; j < n; j += i) {
            is_prime[j] = false;
            meb[j] *= -1;
            if (min_factor[j] == -1) min_factor[j] = i;
            if ((j/i) % i == 0) meb[j] = 0;
        }
    }
    return res;
}

vector<pair<int,int> > prime_factor(int n) {
    vector<pair<int,int> > res;
    while (n != 1) {
        int prime = min_factor[n];
        int exp = 0;
        while (min_factor[n] == prime) {
            ++exp;
            n /= prime;
        }
        res.push_back(make_pair(prime, exp));
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

int main() {
    const long long MOD = 1000000007;
    long long M; cin >> M;
    preprocess();

    long long res = 1;
    for (long long v = 2; v <= M; ++v) {
        if (meb[v] == 0) continue;
        long long r = (M/v) * modinv(M, MOD) % MOD;
        long long add = r * modinv((MOD + 1 - r) % MOD, MOD) % MOD;
        if (meb[v] == -1) res = (res + add) % MOD;
        else res = (res - add + MOD) % MOD;
    }
    cout << res << endl;
}
