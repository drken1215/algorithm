//
// 素因数分解
//
// verified
//   AOJ Course NTL_1_A Elementary Number Theory - Prime Factorize
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=NTL_1_A&lang=jp
//

/*
    入力: n は 2 以上の整数
    出力: 540 = 2^2 × 3^3 × 5 の場合、({2, 2}, {3, 3}, {5, 1})
*/


#include <iostream>
#include <vector>
#include <map>
using namespace std;


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



//------------------------------//
// Examples
//------------------------------//

int main() {
    int n; cin >> n;
    auto res = prime_factorize(n);
    cout << n << ":";
    for (auto prime_beki : res) {
        for (int i = 0; i < prime_beki.second; ++i) {
            cout << " " << prime_beki.first;
        }
    }
    cout << endl;
}



/*
 おまけ
 
 map 形式の出力:
 540 = 2^2 × 3^3 × 5 の場合、{2: 2, 3: 3, 5: 1}
 
map<long long, long long> prime_factorize(long long n) {
    map<long long, long long> res;
    for (long long p = 2; p * p <= n; ++p) {
        while (n % p == 0) { ++res[p]; n /= p; }
    }
    if (n != 1) res[n] = 1;
    return res;
}
*/
