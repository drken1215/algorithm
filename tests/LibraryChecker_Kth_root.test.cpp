#define PROBLEM "https://judge.yosupo.jp/problem/kth_root_integer"

#include "kth_root.hpp"

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    
    int T;
    cin >> T;
    while (T--) {
        uint64_t A, K;
        cin >> A >> K;
        cout << kth_root(A, K) << '\n';
    }
}
