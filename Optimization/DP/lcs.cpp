//
// 最長共通部分列 (LCS) を求める
//
// verifed
//   EDPC F - LCS
//     https://atcoder.jp/contests/dp/tasks/dp_f
//


#include <bits/stdc++.h>
using namespace std;


template<class T> bool chmax(T& a, T b) { if (a < b) { a = b; return true; } return false; }
template<class T> bool chmin(T& a, T b) { if (a > b) { a = b; return true; } return false; }





//------------------------------//
// Examples
//------------------------------//

void EDPC_F() {
    string S, T;
    cin >> S >> T;
    cout << findLCS(S, T) << endl;
}


int main() {
    EDPC_F();
}
