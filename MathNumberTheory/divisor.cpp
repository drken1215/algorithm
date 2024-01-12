//
// 正の整数 n の約数を列挙する, O(√n)
//
// verified
//   ABC 112 D - Partition
//     https://beta.atcoder.jp/contests/abc112/tasks/abc112_d
//


/*
	n の約数を返す
 */


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


vector<long long> calc_divisor(long long n) {
    vector<long long> res;
    for (long long i = 1LL; i*i <= n; ++i) {
        if (n % i == 0) {
            res.push_back(i);
            long long j = n / i;
            if (j != i) res.push_back(j);
        }
    }
    sort(res.begin(), res.end());
    return res;
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    long long N, M;
    cin >> N >> M;
    vector<long long> div = calc_divisor(M);
    
    // M の約数 d であって、d * N <= M となる最大の d を求める
    long long res = 1;
    for (auto d : div) {
        if (d * N <= M) res = max(res, d);
    }
    
    cout << res << endl;
}
