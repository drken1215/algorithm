//
// Monotone 単一始点最短路問題 by D&D Monotone Minima
//   頂点数 N+1 の DAG, 頂点 i, j 間のコスト f(i, j) が Monotone であることを仮定 (argmin が単調非減少)
//   O(N (log N)^2)
//
// verified
//   AtCoder EDPC Z - Frog 3
//     https://atcoder.jp/contests/dp/tasks/dp_z 
//
//   Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
//     https://codeforces.com/contest/319/problem/C 
// 
//   yukicoder No.705 ゴミ拾い Hard
//     https://yukicoder.me/problems/no/705 
//
// Reference:
//   tatyam: Monge の手引き書
//     https://speakerdeck.com/tatyam_prime/monge-noshou-yin-shu
//   


#include <bits/stdc++.h>
using namespace std;


// find min_j f(i, j) for all i, by Monotone Minima, O(H + W log H)
// f(i, j) must be monotone (argmin is not decreasing)
template<class VAL, class FUNC> void MonotoneMinimaRec
(int HL, int HR, int WL, int WR, const FUNC &f, vector<pair<VAL, int>> &res) {
    if (HR - HL <= 0) return;
    int HM = (HL + HR) / 2;
    res[HM].second = WL;
    for (int i = WL; i < WR; i++) {
        VAL val = f(HM, i);
        if (res[HM].first > val) res[HM] = make_pair(val, i);
    }
    MonotoneMinimaRec(HL, HM, WL, res[HM].second + 1, f, res);
    MonotoneMinimaRec(HM + 1, HR, res[HM].second, WR, f, res);
}
template<class VAL, class FUNC> vector<pair<VAL, int>> MonotoneMinima(int H, int W, const FUNC &f) {
    vector<pair<VAL, int>> res(H, make_pair(numeric_limits<VAL>::max() / 2, -1));
    MonotoneMinimaRec(0, H, 0, W, f, res);
    return res;
}

// find shortest path on DAG with monotone cost, by D&D Monotone Minima, O(N (log N)^2)
// vertex: 0, 1, 2, ..., N
// f(i, j) must be monotone (argmin is not decreasing)
template<class VAL, class FUNC> vector<pair<VAL, int>> DDMonotoneMinima(int N, const FUNC &f) {
    vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
    res[0].first = VAL(0);
    auto f2 = [&](int i, int j) -> VAL { return res[j].first + f(j, i); };
    auto rec = [&](auto &&rec, int left, int right) -> void {
        if (right - left <= 1) return;
        int mid = (left + right) / 2;
        rec(rec, left, mid);
        MonotoneMinimaRec(mid, right, left, mid, f2, res);
        rec(rec, mid, right);
    };
    rec(rec, 0, N + 1);
    return res;
}


//------------------------------//
// Examples
//------------------------------//

// AtCoder EDPC Z - Frog 3
/*
    H は単調増加数列
    chmin(dp[j], dp[i] + (H[j] - H[i])^2 + C)
    i -> j のコスト：(H[j] - H[i])^2 ...... 差の凸関数は Monge
    スタート: 0, ゴール: N-1
*/
void EDPC_Z() {
    long long N, C;
    cin >> N >> C;
    vector<long long> H(N);
    for (long long i = 0; i < N; i++) cin >> H[i];
    auto func = [&](int i, int j) -> long long {
        return (H[j] - H[i]) * (H[j] - H[i]) + C;
    };
    auto res = DDMonotoneMinima<long long>(N-1, func);
    cout << res[N-1].first << endl;
}

// Codeforces Round 189 (Div. 1) C. Kalila and Dimna in the Logging Industry
/*
    A: 単調増加, B: 単調減少, ともに長さ N
    i -> j のコストが、B[i] × A[j] で与えられる ..... 単調増加 × 単調減少は Monge
    スタート: 0, ゴール: N-1
*/
void Codeforces_189_C() {
    long long N;
    cin >> N;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> B[i];
    auto func = [&](int i, int j) -> long long {
        return B[i] * A[j];
    };
    auto res = DDMonotoneMinima<long long>(N-1, func);
    cout << res[N-1].first << endl;
}

// yukicoder No.705 ゴミ拾い Hard
/*
    A, X, Y: N 個
    これらを区間に分割していく
    　dp[j] = min_{0 ≦ i < j}(dp[i] + |A[j-1] - X[i]|^3 + |-Y[i]|^3)
    i -> j のコスト：|A[j-1] - X[i]|^3 + |-Y[i]|^3 ...... 差の凸関数 (Monge) + 縞々 (Monge) -> Monge
    スタート: 0, ゴール: N
*/
void yukicoder_705() {
    int N;
    cin >> N;
    vector<long long> A(N), X(N), Y(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> X[i];
    for (int i = 0; i < N; i++) cin >> Y[i];
    auto func = [&](int i, int j) -> long long {
        long long dx = abs(A[j-1] - X[i]), dy = abs(Y[i]);
        return dx * dx * dx + dy * dy * dy;
    };
    auto res = DDMonotoneMinima<long long>(N, func);
    cout << res[N].first << endl;
}


int main() {
    //EDPC_Z();
    //Codeforces_189_C();
    yukicoder_705();
}