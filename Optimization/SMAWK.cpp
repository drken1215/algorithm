//
// Totally Monotone 行最小値問題 by SMAWK
//   H x W 行列の各行の min を求める, Totally Monotone 性を仮定 (列方向に >< がない)
//   O(H + W)
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
// f(i, j) must be totally monotone
template<class VAL, class FUNC> void SMAWKRec
(const vector<int> &X, const vector<int> &Y, const FUNC &f, vector<pair<VAL, int>> &res) {
    if (X.empty()) return;

    // Reduce Step
    vector<int> X2, Y2;
    for (auto y : Y) {
        while (!Y2.empty()) {
            int py = Y2.back(), x = X[(int)Y2.size() - 1];
            if (f(x, y) >= f(x, py)) break;
            Y2.pop_back();
        }
        if (Y2.size() < X.size()) Y2.emplace_back(y);
    }

    // Recurse Step
    for (int i = 1; i < (int)X.size(); i += 2) X2.emplace_back(X[i]);
    SMAWKRec(X2, Y2, f, res);

    // Interpolate Step
    int p = 0;
    for (int i = 0; i < (int)X.size(); i += 2) {
        int lim = (i + 1 < (int)X.size() ? res[X[i + 1]].second : Y.back()), best = Y[p];
        while (Y[p] < lim) {
            p++;
            if (f(X[i], Y[p]) < f(X[i], best)) best = Y[p];
        }
        res[X[i]] = {f(X[i], best), best};
    }
}
template<class VAL, class FUNC> vector<pair<VAL, int>> SMAWK(int H, int W, const FUNC &f) {
    if (H == 0) return {};
    assert(W > 0);
    vector<pair<VAL, int>> res(H, make_pair(numeric_limits<VAL>::max() / 2, -1));
    vector<int> X(H), Y(W);
    for (int i = 0; i < H; i++) X[i] = i;
    for (int j = 0; j < W; j++) Y[j] = j;
    SMAWKRec(X, Y, f, res);
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
    auto res = SMAWK<long long>(N, N, func);
    for (int i = 0; i < N; i++) cout << res[i].first << endl;
}


int main() {
    COLOCON_2018_final_C();
}