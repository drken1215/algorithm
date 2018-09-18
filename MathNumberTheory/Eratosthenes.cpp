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


const int MAX = 10001000;      // 素数判定する最大の数
bool IsPrime[MAX];
vector<int> Era(int n = MAX) {
    vector<int> res;
    IsPrime[0] = false; IsPrime[1] = false;
    for (int i = 2; i < n; ++i) IsPrime[i] = true;
    for (int i = 2; i < n; ++i) {
        if (IsPrime[i]) {
            res.push_back(i);
            for (int j = i*2; j < n; j += i) IsPrime[j] = false;
        }
    }
    return res;
}


int main() {
    vector<int> primes = Era(); // MAX 以下の素数を列挙しておく
    int n;
    while (cin >> n) {
        int num = upper_bound(primes.begin(), primes.end(), n) - primes.begin(); // n 以下が何個か
        cout << num << endl;
    }
}
