//
// GCD (Euclid's algorithm)
//
// cf.
//   拡張ユークリッドの互除法 〜 一次不定方程式 ax + by = c の解き方 〜
//     https://qiita.com/drken/items/b97ff231e43bce50199a
//
// verified
//   ABC 109 C - Skip
//     https://beta.atcoder.jp/contests/abc109/tasks/abc109_c
//


#include <iostream>
#include <vector>
#include <cmath>
using namespace std;


long long GCD(long long a, long long b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    if (b == 0) return a;
    else return GCD(b, a % b);
}



//------------------------------//
// Examples
//------------------------------//

int main() {
    int N; long long X;
    cin >> N >> X;
    vector<long long> x(N);
    for (int i = 0; i < N; ++i) cin >> x[i], x[i] -= X;
    long long g = abs(x[0]);
    for (int i = 0; i < N; ++i) g = GCD(g, abs(x[i]));
    cout << g << endl;
}
