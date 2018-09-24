//
// 重複組合せ nHr = n+r-1Cr
//
// cf.
//   高校数学の美しい物語: 重複組合せの公式と例題（玉，整数解の個数）
//     https://mathtrain.jp/tyohukuc
//
//   よくやる二項係数 (nCk mod. p)、逆元 (a^-1 mod. p) の求め方
//     http://drken1215.hatenablog.com/entry/2018/06/08/210000
//
//
// verified:
//   ARC 110 D - Factorization
//     https://beta.atcoder.jp/contests/abc110/tasks/abc110_d
//


#include <iostream>
#include <vector>
using namespace std;


// 素因数分解
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


// 二項係数
const int MAX = 210000;
const int MOD = 1000000007;

long long fac[MAX], finv[MAX], inv[MAX];
void COMinit(){
    fac[0] = fac[1] = 1;
    finv[0] = finv[1] = 1;
    inv[1] = 1;
    for(int i = 2; i < MAX; i++){
        fac[i] = fac[i-1] * i % MOD;
        inv[i] = MOD - inv[MOD%i] * (MOD/i) % MOD;
        finv[i] = finv[i-1] * inv[i] % MOD;
    }
}
long long COM(int n, int k){
    if(n < k) return 0;
    if (n < 0 || k < 0) return 0;
    return fac[n] * (finv[k] * finv[n-k] % MOD) % MOD;
}


// main
int main() {
    int N, M;
    COMinit();
    cin >> N >> M;
    auto vec = prime_factorize(M);
    long long res = 1;
    for (auto pa : vec) {
        int k = pa.second;
        long long tmp = COM(N+k-1, k);  // nHk = n+k-1Ck
        res = (res * tmp) % MOD;
    }
    cout << res << endl;
}
