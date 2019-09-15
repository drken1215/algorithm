//
// ローリングハッシュ
//
// verified:
//   ABC 141 E - Who Says a Pun?
//     https://atcoder.jp/contests/abc141/tasks/abc141_e
//


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;


// ローリングハッシュ
struct RollingHash {
    static const int base1 = 1007, base2 = 2009;
    static const int mod1 = 1000000007, mod2 = 1000000009;
    vector<long long> hash1, hash2, power1, power2;

    // construct
    RollingHash(const string &S) {
        int n = (int)S.size();
        hash1.assign(n+1, 0);
        hash2.assign(n+1, 0);
        power1.assign(n+1, 1);
        power2.assign(n+1, 1);
        for (int i = 0; i < n; ++i) {
            hash1[i+1] = (hash1[i] * base1 + S[i]) % mod1;
            hash2[i+1] = (hash2[i] * base2 + S[i]) % mod2;
            power1[i+1] = (power1[i] * base1) % mod1;
            power2[i+1] = (power2[i] * base2) % mod2;
        }
    }
    
    // get hash of S[left:right]
    inline pair<long long, long long> get(int l, int r) const {
        long long res1 = hash1[r] - hash1[l] * power1[r-l] % mod1;
        if (res1 < 0) res1 += mod1;
        long long res2 = hash2[r] - hash2[l] * power2[r-l] % mod2;
        if (res2 < 0) res2 += mod2;
        return {res1, res2};
    }

    // get lcp of S[a:] and T[b:]
    inline int getLCP(int a, int b) const {
        int len = min((int)hash1.size()-a, (int)hash1.size()-b);
        int low = 0, high = len;
        while (high - low > 1) {
            int mid = (low + high) >> 1;
            if (get(a, a+mid) != get(b, b+mid)) high = mid;
            else low = mid;
        }
        return low;
    }
};


int main() { 
    int N;
    string S;
    cin >> N >> S;

    // ロリハ
    RollingHash rh(S);

    // 求める
    int res = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
            int lcp = rh.getLCP(i, j);
            lcp = min(lcp, j-i);
            res = max(res, lcp);
        }
    }   
        
    cout << res << endl;
}
