//
// 各辺のコストをスライドアクセスで計算する問題
//
// verified
//   Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem
//     https://codeforces.com/contest/868/problem/F
//
// Reference:
//   noshi: 簡易版 LARSCH Algorithm
//     https://noshi91.hatenablog.com/entry/2023/02/18/005856
//  


#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// よく使うアルゴリズムたち
//------------------------------//

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

// find the d-edges shortest path for d = 1, 2, ..., D, in O(ND)
// vertex: 0, 1, 2, ..., N
// you can choose SMAWK (random access) or Monotone Minima (slide access)
template<class VAL> struct MongeShortestPathWithDEdges {
    VAL INF = numeric_limits<VAL>::max() / 2;

    // solver
    template<class FUNC> vector<VAL> solve(int N, const FUNC &f, int D, const string &solver = "smawk") {
        assert(D <= N);
        vector<VAL> res{VAL(0)}, dp(N + 1, INF);
        dp[0] = 0;
        auto f2 = [&](int i, int j) -> VAL { return (j < i ? dp[j] + f(j, i) : INF); };
        for (int d = 1; d <= D; d++) {
            vector<pair<VAL, int>> tmp;
            if (solver == "smawk") tmp = SMAWK<VAL>(N + 1, N + 1, f2);
            else if (solver == "monotone") tmp = MonotoneMinima<VAL>(N + 1, N + 1, f2);
            for (int i = d; i <= N; i++) dp[i] = tmp[i].first;
            res.emplace_back(dp[N]);
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

// Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem
// by Monge 単一始点 d-辺最短路の d = 1, 2, ..., D における列挙 (by SMAWK 法, in O(ND))
void Codeforces438_F_enum() {
    long long N, K;
    cin >> N >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];

    vector<long long> nums(N + 1, 0);
    long long left = 0, right = 0, score = 0;
    auto add = [&](int i) -> void {
        score += nums[A[i]];
        nums[A[i]]++;
    };
    auto del = [&](int i) -> void {
        nums[A[i]]--;
        score -= nums[A[i]];
    };
    auto f = [&](int i, int j) -> long long {
        while (left < i) del(left++);
        while (i < left) add(--left);
        while (right < j) add(right++);
        while (j < right) del(--right);
        return score;
    };
    MongeShortestPathWithDEdges<long long> msp;
    auto res = msp.solve(N, f, K, "monotone");  // we should use slide access
    cout << res[K] << endl;
}


int main() {
    Codeforces438_F_enum();
}