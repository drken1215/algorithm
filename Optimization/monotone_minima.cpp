//
// Monotone 行最小値問題 by Monotone Minima
//   H x W 行列の各行の min を求める, Monotone 性を仮定 (argmin が単調非減少)
//   O(H + W log H)
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//
// Reference:
//   tatyam: Monge の手引き書
//     https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu
//   


#include <bits/stdc++.h>
using namespace std;


// find min_j f(i, j) for all i, by Monotone Minima, O(H + W log H)
// f(i, j) must be monotone (argmin is not decreasing)
template<class VAL, class FUNC> vector<pair<VAL, int>> MonotoneMinima(int H, int W, const FUNC &f) {
    vector<pair<VAL, int>> res(H, make_pair(numeric_limits<VAL>::max() / 2, -1));
    auto rec = [&](auto &&rec, int HL, int HR, int WL, int WR) -> void {
        if (HR - HL <= 0) return;
        int HM = (HL + HR) / 2;
        res[HM].second = WL;
        for (int i = WL; i < WR; i++) {
            VAL val = f(HM, i);
            if (res[HM].first > val) res[HM] = make_pair(val, i);
        }
        rec(rec, HL, HM, WL, res[HM].second + 1);
        rec(rec, HM + 1, HR, res[HM].second, WR);
    };
    rec(rec, 0, H, 0, W);
    return res;
}


//------------------------------//
// Examples
//------------------------------//

// COLOCON 2018 Final C - スペースエクスプローラー高橋君
/*
    M[i][j] = (j - i)^2 は Monge
    列 j に一様に A[j] を加えても Monge
*/
void COLOCON_2018_final_C() {
    long long N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    auto func = [&](long long i, long long j) -> long long {
        return A[j] + (j - i) * (j - i);
    };
    auto res = MonotoneMinima<long long>(N, N, func);
    for (int i = 0; i < N; i++) cout << res[i].first << endl;
}


int main() {
    COLOCON_2018_final_C();
}