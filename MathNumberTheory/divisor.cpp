//
// ³‚Ì®” n ‚Ì–ñ”‚ğ—ñ‹“‚·‚é, O(ãn)
//
// verified
//   ABC 112 D - Partition
//     https://beta.atcoder.jp/contests/abc112/tasks/abc112_d
//


/*
	n ‚Ì–ñ”‚ğ•Ô‚·
*/


#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


vector<long long> calc_divisor(long long n) {
	vector<long long> res;
	for (long long i = 1LL; i*i <= n; ++i) {
		if (n % i == 0) {
			res.push_back(i);
			long long j = n / i;
			if (j != i) res.push_back(j);
		}
	}
	sort(res.begin(), res.end());
	return res;
}


int main() {
	long long N, M;
	cin >> N >> M;
	vector<long long> div = calc_divisor(M);
	
	// M ‚Ì–ñ” d ‚Å‚ ‚Á‚ÄAd * N <= M ‚Æ‚È‚éÅ‘å‚Ì d ‚ğ‹‚ß‚é
	long long res = 1;
	for (auto d : div) {
		if (d * N <= M) res = max(res, d);
	}

	cout << res << endl;
}