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


struct Manacher {
    string S;
    vector<int> len;

    // construct
    Manacher(const string &iS) { init(iS); }
    void init(const string &iS) {
        S = "";
        for (int i = 0; i < (int)iS.size(); ++i) {
            S += iS[i];
            if (i + 1 < (int)iS.size()) S += "$";
        }
        construct();
    }
    void construct() {
        int N = (int)S.size();
        len.resize(N);
        int i = 0, j = 0;
        while (i < N) {
            while (i-j >= 0 && i+j < N && S[i-j] == S[i+j]) ++j;
            len[i] = j;
            int k = 1;
            while (i-k >= 0 && i+k < N && k+len[i-k] < j) len[i+k] = len[i-k], ++k;
            i += k, j -= k;
        }
    }

    // radius, center is i
    inline int get_odd(int i) { return (len[i*2] + 1) / 2; }

    // radius, center is between i-1 and i
    inline int get_even(int i) { return len[i*2-1] / 2; }

    // judge if [left, right) is palindrome
    inline bool ispalin(int left, int right) {
        int mid = (left + right) / 2;
        if ((right - left) & 1) return ( get_odd(mid) == (right - left + 1)/2);
        else return (get_even(mid) == (right - left)/2);
    }
};


int main() {
    int N; string S;
    cin >> N >> S;

    Manacher m(S);
    int res = N;
    for (int len = 2; len < N; ++len) {
        bool ok = true;
        for (int i = len - 1; i < N; i += len - 1) {
            if (m.get_odd(i) != min(i+1, N-i)) ok = false;
        }
        if (ok) {
            res = len;
            break;
        }
    }
    cout << res << endl;
}
