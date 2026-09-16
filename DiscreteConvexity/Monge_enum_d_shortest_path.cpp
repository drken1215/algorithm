//
// Monge 単一始点 d-辺最短路の d = 1, 2, ..., D における列挙 (by SMAWK 法, in O(ND))
//
// verified
//   yukicoder No.952 危険な火薬庫
//     https://yukicoder.me/problems/no/952
//
//   Monoxer Programming Contest for Engineers K - Coupon
//     https://mofecoder.com/contests/monoxercon_202508/tasks/monoxercon_202508_k
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

// find the d-edges shortest path for d = 1, 2, ..., D, in O(ND)
// vertex: 0, 1, 2, ..., N
// you can choose SMAWK or Monotone Minima
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
            for (int i = d; i <= N; i++) dp[i] = tmp[i].first;
            res.emplace_back(dp[N]);
        }
        return res;
    }
};


//------------------------------//
// Examples
//------------------------------//

// yukicoder No.952 危険な火薬庫
/*
    末尾にドアを追加する。
    ノード 0, 1, ..., N, N+1
    f(i, j) = (S[j-1] - S[i])^2
    選ぶドアが k 個 → 区間の個数は、N+1-k
*/
void yukicoder_952() {
    int N;
    cin >> N;
    vector<long long> A(N + 1, 0), S(N + 2, 0);
    for (int i = 0; i < N; i++) cin >> A[i], S[i + 1] = S[i] + A[i];
    auto f = [&](int i, int j) -> long long { 
        return (S[j - 1] - S[i]) * (S[j - 1] - S[i]);
    };
    MongeShortestPathWithDEdges<long long> msp;
    auto res = msp.solve(N + 1, f, N + 1);
    for (int k = 1; k <= N; k++) cout << res[N + 1 - k] << endl;
}

// Monoxer Programming Contest for Engineers K - Coupon
/*
    f(i, j) = (S[j] - S[i]) - max_{i ≦ k < j} A[k] / 2 + (j - i - 1)^2
*/
template<class MeetSemiLattice> struct SparseTable {
    using Func = function<MeetSemiLattice(MeetSemiLattice, MeetSemiLattice)>;

    // core member
    Func OP = [](const MeetSemiLattice &l, const MeetSemiLattice &r) {
        return min(l, r);
    };
    vector<vector<MeetSemiLattice>> dat;
    vector<int> height;
    
    SparseTable() {}
    SparseTable(const vector<MeetSemiLattice> &vec) {
        init(vec);
    }
    SparseTable(const vector<MeetSemiLattice> &vec, const Func &op)  {
        init(vec, op);
    }
    void init(const vector<MeetSemiLattice> &vec) {
        int n = (int)vec.size(), h = 1;
        while ((1<<h) <= n) ++h;
        dat.assign(h, vector<MeetSemiLattice>(1<<h));
        height.assign(n+1, 0);
        for (int i = 2; i <= n; i++) height[i] = height[i>>1]+1;
        for (int i = 0; i < n; ++i) dat[0][i] = vec[i];
        for (int i = 1; i < h; ++i) {
            for (int j = 0; j < n; ++j)
                dat[i][j] = OP(dat[i-1][j], dat[i-1][min(j+(1<<(i-1)),n-1)]);
        }
    }
    void init(const vector<MeetSemiLattice> &vec, const Func &op) {
        OP = op;
        init(vec);
    }
    
    MeetSemiLattice get(int a, int b) {
        return OP(dat[height[b-a]][a], dat[height[b-a]][b-(1<<height[b-a])]);
    }
};
void MonoxerCoupon() {
    long long N;
    cin >> N;
    vector<long long> A(N), S(N+1, 0);
    for (int i = 0; i < N; i++) cin >> A[i], S[i+1] = S[i] + A[i];
    SparseTable<long long> sparsetable(A, [](long long p, long long q){return max(p, q);});
    auto f = [&](int i, int j) -> long long {
        return S[j] - S[i] - sparsetable.get(i, j)/2 + (j - i - 1) * (j - i - 1);
    };
    MongeShortestPathWithDEdges<long long> msp;
    auto res = msp.solve(N, f, N);
    for (int i = 1; i <= N; i++) cout << res[i] << endl;
}


int main() {
    //yukicoder_952();
    MonoxerCoupon();
}