#define PROBLEM "https://atcoder.jp/contests/abc361/tasks/abc361_f"

#include "kth_root.hpp"

int main() {
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
