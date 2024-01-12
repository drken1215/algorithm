//
// nCk mod m、オーソドックス
//
// cf.
//   drken: 「1000000007 で割ったあまり」の求め方を総特集！ 〜 逆元から離散対数まで 〜
//     https://qiita.com/drken/items/3b4fdf0a78e7a138cd9a
//
// verified:
//   AGC 025 B - RGB Coloring
//     https://beta.atcoder.jp/contests/agc025/tasks/agc025_b
//


#include <iostream>
using namespace std;


const int MAX = 510000;
const int MOD = 998244353;

long long fac[MAX], finv[MAX], inv[MAX];
void COMinit(){
    fac[0] = fac[1] = 1;
    finv[0] = finv[1] = 1;
    inv[1] = 1;
    for (int i = 2; i < MAX; i++){
        fac[i] = fac[i - 1] * i % MOD;
        inv[i] = MOD - inv[MOD%i] * (MOD / i) % MOD;
        finv[i] = finv[i - 1] * inv[i] % MOD;
    }
}
long long COM(int n, int k){
    if (n < k) return 0;
    if (n < 0 || k < 0) return 0;
    return fac[n] * (finv[k] * finv[n - k] % MOD) % MOD;
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    COMinit(); // テーブルを作る
    long long N, A, B, K;
    cin >> N >> A >> B >> K;
    long long res = 0;
    for (long long a = 0; a <= N; ++a) { // A の個数を全探索
        long long rem = K - a * A;
        if (rem % B != 0) continue;
        long long b = rem / B;
        if (b > N) continue;
        long long tmp = COM(N, a) * COM(N, b) % MOD;
        res += tmp;
        res %= MOD;
    }
    
    cout << res << endl;
}
