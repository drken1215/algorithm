//
// nCk mod m、漸化式
//
// cf.
//   よくやる二項係数 (nCk mod. p)、逆元 (a^-1 mod. p) の求め方
//     http://drken1215.hatenablog.com/entry/2018/06/08/210000
//
// verified:
//   CODECHEF July Challenge 2018 Division 2 D - No Minimum No Maximum
//     https://www.codechef.com/problems/NMNMX
//


/*
    nCk mod m
 
    ・n も k も 5000 程度と小さいとき、動的計画法によるテーブル生成が有力
    ・m が素数でなくてもよいのが大きな魅力 (下の例では m = 1000000006)
*/


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


// 二項係数
const int MOD = 1000000006;
const int MAX_C = 5100;
long long Com[MAX_C][MAX_C];
void COMinit() {
    Com[0][0] = 1;
    for (int i = 1; i < MAX_C; ++i) Com[0][i] = 0;
    for (int i = 1; i < MAX_C; ++i) {
        Com[i][0] = 1;
        for (int j = 1; j < MAX_C; ++j) {
            Com[i][j] = (Com[i-1][j-1] + Com[i-1][j]) % MOD;
        }
    }
}
long long COM(int n, int k){
    if (n < k) return 0;
    if (n < 0 || k < 0) return 0;
    return Com[n][k];
}


// mod 累乗
long long modpow(long long a, long long n, long long mod) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}


int main() {
    COMinit();
    const int MOD2 = 1000000007;
    int T; cin >> T;
    for (int _ = 0; _ < T; ++_) {
        long long N, K; cin >> N >> K;
        vector<long long> a(N);
        for (int i = 0; i < N; ++i) {
            cin >> a[i];
        }
        sort(a.begin(), a.end());
        
        long long res = 1;
        for (long long i = 1; i + 1 < N; ++i) {
            long long left = COM(N-i-1, K-1);
            long long right = COM(i, K-1);
            long long all = COM(N-1, K-1);
            long long count = (all - left - right + MOD*2) % MOD;
            long long fact = modpow(a[i], count, MOD2);
            res = res * fact % MOD2;
        }
        cout << res << endl;
    }
}
