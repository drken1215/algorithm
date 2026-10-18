//
// Monge 最短路問題に関するアルゴリズム全集
//
// verified
//   ABC 480 G - Zonal Score Minimization
//     https://atcoder.jp/contests/abc480/tasks/abc480_g
//


#include <bits/stdc++.h>
using namespace std;


// find min_j f(i, j) for all i, by , O(H + W log H)
// f(i, j) must be totally monotone
// slide access NG
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

// find shortest path on DAG with anti-Monge cost, by D&D SMAWK, O(NK log N)
// vertex: 0, 1, 2, ..., N
// f(i, j) must be anti-Monge, j >= i + K
template<class VAL, class FUNC> vector<pair<VAL, int>> antiMongeDDSMAWK(int N, const FUNC &f, int K = 1) {
    assert(K >= 1);
    const VAL INF = numeric_limits<VAL>::max() / 4;
    vector<pair<VAL, int>> res(N + 1, make_pair(INF, -1));
    vector<pair<VAL, int>> tmp(N + 1, {INF, -1});
    res[0].first = VAL(0);
    auto relax = [&](int i, int j) -> void {
        if (i != 0 && res[i].second == -1) return;
        VAL val = res[i].first + f(i, j);
        if (val < res[j].first) res[j] = {val, i};
    };
    auto rec = [&](auto &&rec, int left, int right) -> void {
        if (right - left <= 1) return;
        int mid = (left + right) / 2, border = min(right, mid + K - 1);
        rec(rec, left, mid);
        for (int j = mid; j < border; j++) {
            for (int i = left; i < mid && i <= j - K; i++) relax(i, j);
        }
        if (border < right) {
            vector<int> X, Y;
            for (int j = border; j < right; j++) X.push_back(j);
            for (int i = mid - 1; i >= left; i--) {
                if (i == 0 || res[i].second != -1) Y.push_back(i);
            }
            if (!Y.empty()) {
                vector<int> cols(Y.size());
                iota(cols.begin(), cols.end(), 0);
                auto f2 = [&](int j, int k) -> VAL { return res[Y[k]].first + f(Y[k], j); };
                SMAWKRec(X, cols, f2, tmp);
                for (auto x : X) {
                    if (tmp[x].first < res[x].first) res[x] = {tmp[x].first, Y[tmp[x].second]};
                }
            }
        }
        rec(rec, mid, right);
    };
    rec(rec, 0, N + 1);
    return res;
}


//------------------------------//
// Examples
//------------------------------//

// ABC 480 G - Zonal Score Minimization
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
void ABC_480_G() {
    int N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    SparseTable<long long> stmi(A);
    SparseTable<long long> stma(A, [&](long long p, long long q){return max(p, q);});

    auto func = [&](int l, int r) -> long long {
        assert(r - l >= 2);
        return stma.get(l, r) - stmi.get(l, r);
    };
    auto res = antiMongeDDSMAWK<long long>(N, func, 2);
    cout << res.back().first << endl;
}

                  
int main() {
    ABC_480_G();
}