//
// 1：Monotone 単一始点最短路問題 by D&D Monotone Minima
//   頂点数 N+1 の DAG, 頂点 i, j 間のコスト f(i, j) が Monotone であることを仮定 (argmin が単調非減少)
//   O(N (log N)^2)
//
// 2：行最小値問題 by Monotone Minima
//   H x W 行列の各行の min を求める, Monotone 性を仮定 (argmin が単調非減少)
//   O(H + W log H)
//
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//
//   AtCoder EDPC Z - Frog 3
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//
//   Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
//     https://codeforces.com/contest/319/problem/C 
//
//
// Reference:
//   tatyam: Monge の手引き書
//     https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu
//   
//   noshi: 簡易版 LARSCH Algorithm
//     https://noshi91.hatenablog.com/entry/2023/02/18/005856
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

// find shortest path on DAG with monotone cost, by D&D Monotone Minima, O(N (log N)^2)
// vertex: 0, 1, 2, ..., N
// f(i, j) must be monotone (argmin is not decreasing)
template<class VAL, class FUNC> vector<pair<VAL, int>> MonotoneMinimaDD(int N, const FUNC &f) {
    vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
    res[0].first = VAL(0);
    auto f2 = [&](int i, int j) -> VAL { return res[j].first + f(j, i); };
    auto rec2 = [&](auto &&rec2, int HL, int HR, int WL, int WR) -> void {
        if (HR - HL <= 0) return;
        int HM = (HL + HR) / 2;
        res[HM].second = WL;
        for (int i = WL; i < WR; i++) {
            VAL val = f2(HM, i);
            if (res[HM].first > val) res[HM] = make_pair(val, i);
        }
        rec2(rec2, HL, HM, WL, res[HM].second + 1);
        rec2(rec2, HM + 1, HR, res[HM].second, WR);
    };
    auto rec1 = [&](auto &&rec1, int left, int right) -> void {
        if (right - left <= 1) return;
        int mid = (left + right) / 2;
        rec1(rec1, left, mid);
        rec2(rec2, mid, right, left, mid);
        rec1(rec1, mid, right);
    };
    rec1(rec1, 0, N + 1);
    return res;
}


//------------------------------//
// Examples
//------------------------------//

// COLOCON 2018 Final C - スペースエクスプローラー高橋君
void COLOCON_2018_final_C() {
    long long N;
    cin >> N;
    vector<long long> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    auto func = [&](int i, int j) -> long long {
        return (long long)(i - j) * (i - j) + a[j];
    };
    auto res = MonotoneMinima<long long>(N, N, func);
    for (int i = 0; i < N; ++i) cout << res[i].first << endl;
}

// AtCoder EDPC Z - Frog 3
void EDPC_Z() {
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N);
    for (long long i = 0; i < N; i++) cin >> H[i];
    auto func = [&H, &C](int i, int j) -> long long {
        return (H[j] - H[i]) * (H[j] - H[i]) + C;
    };
    auto res = MonotoneMinimaDD<long long>(N-1, func);
    cout << res[N-1].first << endl;
}

// Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
void Codeforces_189_C() {
    int N;
    cin >> N;
    vector<long long> a(N), b(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    for (int i = 0; i < N; ++i) cin >> b[i];
    auto func = [&](int i, int j) -> long long { return a[j] * b[i]; };
    auto dp = MonotoneMinimaDD<long long>(N-1, func);
    cout << dp[N-1].first << endl;
}


int main() {
    COLOCON_2018_final_C();
    //EDPC_Z();
    //Codeforces_189_C();
}