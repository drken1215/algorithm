//
// ローリングハッシュ
//
// verified:
//   ARC 060 F - 最良表現
//     https://atcoder.jp/contests/arc060/tasks/arc060_d
//


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;


struct RollingHash {
    const int base = 9973;
    const vector<int> mod = {999999937LL, 1000000007LL};
    string S_;
    vector<long long> hash[2], power[2];
    RollingHash(){}
    RollingHash(const string &S) : S_(S) {
        int n = (int)S.size();
        for (int iter = 0; iter < 2; ++iter) {
            hash[iter].assign(n+1, 0);
            power[iter].assign(n+1, 1);
            for (int i = 0; i < n; ++i) {
                hash[iter][i+1] = (hash[iter][i] * base + S[i]) % mod[iter];
                power[iter][i+1] = power[iter][i] * base % mod[iter];
            }
        }
    }
    // get hash of S[left:right]
    inline long long get(int l, int r, int id = 0) const {
        long long res = hash[id][r] - hash[id][l] * power[id][r-l] % mod[id];
        if (res < 0) res += mod[id];
        return res;
    }
    // get lcp of S[a:] and S[b:]
    inline int getLCP(int a, int b) const {
        int len = min((int)S_.size()-a, (int)S_.size()-b);
        int low = -1, high = len + 1;
        while (high - low > 1) {
            int mid = (low + high) / 2;
            if (get(a, a+mid, 0) != get(b, b+mid, 0)) high = mid;
            else if (get(a, a+mid, 1) != get(b, b+mid, 1)) high = mid;
            else low = mid;
        }
        return low;
    }
    // get lcp of S[a:] and T[b:]
    inline int getLCP(const RollingHash &t, int a, int b) const {
        int len = min((int)S_.size()-a, (int)S_.size()-b);
        int low = -1, high = len + 1;
        while (high - low > 1) {
            int mid = (low + high) / 2;
            if (get(a, a+mid, 0) != get(b, b+mid, 0)) high = mid;
            else if (get(a, a+mid, 1) != get(b, b+mid, 1)) high = mid;
            else low = mid;
        }
        return low;
    }
};


vector<long long> divisor(long long n) {
    vector<long long> res;
    for (long long i = 1LL; i*i <= n; ++i) {
        if (n%i == 0LL) {
            res.push_back(i);
            long long temp = n/i;
            if (i != temp) res.push_back(temp);
        }
    }
    sort(res.begin(), res.end());
    return res;
}
    
int main() {
    string str; cin >> str;
    int n = (int)str.size();
    vector<long long> divs = divisor(n);
    long long syuuki = n;
    for (auto d : divs) {
        bool ok = true;
        for (int j = 0; j + d < n; ++j) {
            if (str[j] != str[j+d]) ok = false;
        }
        if (ok) syuuki = min(syuuki, d);
    }
    if (syuuki == n) cout << 1 << endl << 1 << endl;
    else if (syuuki == 1) cout << n << endl << 1 << endl;
    else {
        vector<int> cannot_cut(n*2, 0);
        RollingHash rh(str);
        for (int d = 1; d < n; ++d) {
            if (cannot_cut[d]) continue;
            for (int dd = d; dd < n; dd += d) {
                if (rh.get(0, d) != rh.get(dd, dd+d)) break;
                cannot_cut[dd + d] = true;
            }
            for (int dd = n-d*2; dd >= 0; dd -= d) {
                if (rh.get(dd, dd+d) != rh.get(n-d, n)) break;
                cannot_cut[dd] = true;
            }
        }
        int con = 0;
        for (int i = 1; i < n; ++i) if (!cannot_cut[i]) ++con;
        cout << 2 << endl << con << endl;
    }
}
