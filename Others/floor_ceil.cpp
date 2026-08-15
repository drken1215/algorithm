//
// 負数にも対応した floor, ceil
//
// verified:
//   ABC 334 B - Christmas Trees
//     https://atcoder.jp/contests/abc334/tasks/abc334_b
//


#include <bits/stdc++.h>
using namespace std;


// floor, ceil
template<class T> T floor(T a, T b) {
    if (a % b == 0 || a >= 0) return a / b;
    else return -((-a) / b) - 1;
}
template<class T> T ceil(T x, T y) {
    return floor(x + y - 1, y);
}


//------------------------------//
// Examples
//------------------------------//

// ABC 334 B - Christmas Trees
void ABC_334_B() {
    long long A, M, L, R;
    cin >> A >> M >> L >> R;
    long long l = ceil(L - A, M), r = floor(R - A, M);
    cout << r - l + 1 << endl;
}


int main() {
    ABC_334_B();
}