//
// 高速ゼータ変換 (と高速メビウス変換)
//
// verified:
//   AtCoder ARC 100 E - Or Plus Max
//     https://atcoder.jp/contests/arc100/tasks/arc100_c
//     解説：https://drken1215.hatenablog.com/entry/2018/07/25/200000
//
//


// ARC 100 E の解答例
#include <iostream>
#include <vector>
using namespace std;

typedef pair<long long, long long> pll;
const long long INF = 1LL<<60;

void chmax(pll &a, pll b) {
  if (a.first < b.first) a.second = max(a.first, b.second), a.first = b.first;
  else a.second = max(a.second, b.first);
}

int main() {
  int N; cin >> N;
  vector<long long> A(1<<N); vector<pll> dp(1<<N);
  for (int bit = 0; bit < (1<<N); ++bit) {
    cin >> A[bit];
    dp[bit] = {A[bit], -INF};
  }
    
  // 高速ゼータ変換
  for (int i = 0; i < N; ++i)
    for (int bit = 0; bit < (1<<N); ++bit)
      if (bit & (1<<i))
        chmax(dp[bit], dp[bit^(1<<i)]);

  // 答えを求める
  long long res = -INF;
  for (int bit = 0; bit < (1<<N); ++bit) {
    long long tmp = dp[bit].first + dp[bit].second;
    res = max(res, tmp);
    if (bit) cout << res << endl;
  }
}
