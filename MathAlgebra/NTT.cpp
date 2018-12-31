//
// NTT (Number Theoretic Transform)
//
// cf.
//
//
// verified:
//   Educational Codeforces Round 057 G - Lucky Tickets
//     https://codeforces.com/contest/1096/problem/G
//

#include <iostream>
#include <vector>
using namespace std;



long long modpow(long long a, long long n, long long mod) {
    long long res = 1;
    while (n > 0) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}

long long modinv(long long a, long long mod) {
    long long b = mod, u = 1, v = 0;
    while (b) {
        long long t = a/b;
        a -= t*b; swap(a, b);
        u -= t*v; swap(u, v);
    }
    u %= mod;
    if (u < 0) u += mod;
    return u;
}

namespace NTT {
    const int MOD = 998244353;  // to be set appropriately
    const long long PR = 3;     // to be set appropriately
    
    void trans(vector<long long> &v, bool inv = false) {
        int n = (int)v.size();
        for (int i = 0, j = 1; j < n-1; j++) {
            for (int k = n>>1; k > (i ^= k); k >>= 1);
            if (i > j) swap(v[i], v[j]);
        }
        for (int t = 2; t <= n; t <<= 1) {
            long long bw = modpow(PR, (MOD-1)/t, MOD);
            if (inv) bw = modinv(bw, MOD);
            for (int i = 0; i < n; i += t) {
                long long w = 1;
                for (int j = 0; j < t/2; ++j) {
                    int j1 = i + j, j2 = i + j + t/2;
                    long long c1 = v[j1], c2 = v[j2] * w % MOD;
                    v[j1] = c1 + c2;
                    v[j2] = c1 - c2 + MOD;
                    while (v[j1] >= MOD) v[j1] -= MOD;
                    while (v[j2] >= MOD) v[j2] -= MOD;
                    w = w * bw % MOD;
                }
            }
        }
        if (inv) {
            long long inv_n = modinv(n, MOD);
            for (int i = 0; i < n; ++i) v[i] = v[i] * inv_n % MOD;
        }
    }
    
    // C is A*B
    vector<long long> mult(vector<long long> A, vector<long long> B) {
        int size_a = 1; while (size_a < A.size()) size_a <<= 1;
        int size_b = 1; while (size_b < B.size()) size_b <<= 1;
        int size_fft = max(size_a, size_b) << 1;
        
        vector<long long> cA(size_fft, 0), cB(size_fft, 0), cC(size_fft, 0);
        for (int i = 0; i < A.size(); ++i) cA[i] = A[i];
        for (int i = 0; i < B.size(); ++i) cB[i] = B[i];
        
        trans(cA); trans(cB);
        for (int i = 0; i < size_fft; ++i) cC[i] = cA[i] * cB[i] % MOD;
        trans(cC, true);
        
        vector<long long> res((int)A.size() + (int)B.size() - 1);
        for (int i = 0; i < res.size(); ++i) res[i] = cC[i];
        return res;
    }
};



int main() {
    int n, k; cin >> n >> k;
    n /= 2;
    vector<long long> v(1<<22, 0);
    for (int i = 0; i < k; ++i) {
        int d; cin >> d;
        v[d] = 1;
    }
    NTT::trans(v);
    vector<long long> each(1<<22, 1);
    for (int i = 0; i < each.size(); ++i) {
        int N = n;
        while (N > 0) {
            if (N & 1) each[i] = each[i] * v[i] % NTT::MOD;
            v[i] = v[i] * v[i] % NTT::MOD;
            N >>= 1;
        }
    }
    NTT::trans(each, true);
    
    long long res = 0;
    for (int i = 0; i < each.size(); ++i) {
        res += each[i] * each[i] % NTT::MOD;
    }
    res %= NTT::MOD;
    cout << res << endl;
}
