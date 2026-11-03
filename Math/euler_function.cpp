//
// Euler 関数
//
// verified
//   AOJ Course NTL_1_D: Euler's Phi Function
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=NTL_1_D&lang=ja
//


#include <iostream>
#include <vector>
using namespace std;

long long Euler(long long N) {
    long long res = N;
    for (long long p = 2; p * p <= N; ++p) {
        if (N % p == 0) {
            res = res / p * (p - 1);
            while (N % p == 0) N /= p;
        }
    }
    if (N != 1) res = res / N * (N - 1);
    return res;
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    long long N;
    cin >> N;
    cout << Euler(N) << endl;
}
