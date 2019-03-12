//
// Manacher のアルゴリズム
//   res[i] := S[i] を中心とする S 上の奇数長の最長回文の半径 (中心を含む)
//   偶数長の回文を求めたいときは S = "abc" に対して "a$b$c" といった処理をするテクがある
//
// verified:
//   RUPC 2019 day3-E 往復文字列
//     https://onlinejudge.u-aizu.ac.jp/beta/room.html#RitsCamp19Day3/problems/E
//


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;


vector<int> Manacher(const string &S) {
    int N = (int)S.size();
    vector<int> res(N);
    int i = 0, j = 0;
    while (i < N) {
        while (i-j >= 0 && i+j < N && S[i-j] == S[i+j]) ++j;
        res[i] = j;
        int k = 1;
        while (i-k >= 0 && i+k < N && k+res[i-k] < j) res[i+k] = res[i-k], ++k;
        i += k, j -= k;
    }
    return res;
}


int main() {
    int N; string S;
    cin >> N >> S;

    auto rad = Manacher(S);
    int res = N;
    for (int len = 2; len < N; ++len) {
        bool ok = true;
        for (int i = len - 1; i < N; i += len - 1) {
            if (rad[i] != min(i+1, N-i)) ok = false;
        }
        if (ok) {
            res = len;
            break;
        }
    }
    cout << res << endl;
}
