//
// 二分探索
//
// reference:
//   drken: 二分探索アルゴリズムを一般化 〜 めぐる式二分探索法のススメ 〜
//　　　https://qiita.com/drken/items/97e37dd6143e33a64c8c
//
// verified
//   競プロ典型 90 問 001 - Yokan Party（★4）
//     https://atcoder.jp/contests/typical90/tasks/typical90_a
//


#include <bits/stdc++.h>
using namespace std;


// binary search
template <class T, class FUNC> T binary_search_01(T ng, T ok, const FUNC &f) {
    while (abs(ok - ng) > 1) {
        T mid = ok + ng >> 1;
        if (f(mid)) ok = mid;
        else ng = mid;
    }
    return ok;
}
template <class T, class FUNC> T binary_search_10(T ok, T ng, const FUNC &f) {
    return binary_search_01(ng, ok, f);
}


//------------------------------//
// Examaple
//------------------------------//

// 競プロ典型 90 問 001 - Yokan Party（★4）
void KyoproTenkei90_001() {
    long long N, L, K;
    cin >> N >> L >> K;
    vector<long long> A(N+2, 0);
    for (int i = 0; i < N; i++) cin >> A[i+1];
    A[N+1] = L;

    // 長さ x 以上ずつ切って、K+1 個を切り出せるか
    auto check = [&](long long x) -> bool {
        int prev = 0;
        long long num = 0;
        for (int i = 0; i < N+2; i++) {
            if (A[i] - A[prev] >= x) {
                prev = i;
                num++;
            }
        }
        return num >= K+1;
    };
    long long res = binary_search_10(1LL, L, check);
    cout << res << endl;
}


int main() {
    KyoproTenkei90_001();
}