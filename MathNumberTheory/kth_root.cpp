//
// 正の整数 N, K に対して、x^K <= N を満たす最大の正の整数 x を求める
//
// verified:
//   Yosupo Library Checker - Kth Root (Integer)
//     https://judge.yosupo.jp/problem/kth_root_integer
//
//   AtCoder ABC 361 F - x = a^b
//     https://atcoder.jp/contests/abc361/tasks/abc361_f
//


#include <bits/stdc++.h>
using namespace std;


// N < 2^64, K <= 64
uint64_t kth_root(uint64_t N, uint64_t K) {
    assert(K >= 1);
    if (N <= 1 || K == 1) return N;
    if (K >= 64) return 1;
    if (N == uint64_t(-1)) --N;
    
    auto mul = [&](uint64_t x, uint64_t y) -> uint64_t {
        if (x < UINT_MAX && y < UINT_MAX) return x * y;
        if (x == uint64_t(-1) || y == uint64_t(-1)) return uint64_t(-1);
        return (x <= uint64_t(-1) / y ? x * y : uint64_t(-1));
    };
    auto power = [&](uint64_t x, uint64_t k) -> uint64_t {
        if (k == 0) return 1ULL;
        uint64_t res = 1ULL;
        while (k) {
            if (k & 1) res = mul(res, x);
            x = mul(x, x);
            k >>= 1;
        }
        return res;
    };
    
    uint64_t res;
    if (K == 2) res = sqrtl(N) - 1;
    else if (K == 3) res = cbrt(N) - 1;
    else res = pow(N, nextafter(1 / double(K), 0));
    while (power(res + 1, K) <= N) ++res;
    return res;
}



//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Kth Root (Integer)
void Yosupo_Kth_Root() {
    int T;
    cin >> T;
    while (T--) {
        uint64_t A, K;
        cin >> A >> K;
        cout << kth_root(A, K) << '\n';
    }
}

// AtCoder ABC 361 F - x = a^b
void ABC_361_F() {
    vector<long long> prs = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59};
    long long N;
    cin >> N;
    
    long long res = 0;
    for (int bit = 1; bit < (1 << prs.size()); ++bit) {
        long long K = 1;
        for (int i = 0; i < prs.size(); ++i) {
            if (bit & (1 << i)) {
                K *= prs[i];
                if (K > 60) break;
            }
        }
        if (__builtin_popcount(bit) % 2 == 1) res += kth_root(N, K);
        else res -= kth_root(N, K);
    }
    cout << res << '\n';
}
    

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    
    //Yosupo_Kth_Root();
    ABC_361_F();
}
