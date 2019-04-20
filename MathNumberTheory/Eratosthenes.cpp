//
// エラトステネスの篩
//
// cf.
//   高校数学の美しい物語: エラトスネテスの篩
//     https://mathtrain.jp/eratosthenes
//
// verified
//   AOJ 0009 素数
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=0009&lang=jp
//

/*
 
    n 以下の素数をすべて列挙する
 
*/


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


vector<bool> isprime;
vector<int> Era(int n) {
    isprime.resize(n, true);
    vector<int> res;
    isprime[0] = false; isprime[1] = false;
    for (int i = 2; i < n; ++i) isprime[i] = true;
    for (int i = 2; i < n; ++i) {
        if (isprime[i]) {
            res.push_back(i);
            for (int j = i*2; j < n; j += i) isprime[j] = false;
        }
    }
    return res;
}


int main() {
    vector<int> primes = Era(1000000);
    int n;
    while (cin >> n) {
        int num = upper_bound(primes.begin(), primes.end(), n) - primes.begin(); // n 以下が何個か
        cout << num << endl;
    }
}
