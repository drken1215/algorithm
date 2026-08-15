//
// Alien's Trick --- min: f(x) s.t. g(x) = K (by ラグランジュ緩和)
// 
/*
    min_{x}: f(x) s.t. g(x) = K
    = max_{λ は整数}: min_{x} (f(x) + λ(g(x) - K))

    【仮定】
　　・f(x), g(x): 0 以上の整数値をとる関数
　　・K: 整数
　　・g(x) = p となる x についての f(x) の最小値を h(p) とおくと、h(p) は p について下に凸な関数
　　・h(K) < ∞ (g(x) = K となる x が存在)

　　 【コメント】
　　・制約条件が g(x) <= K の場合、h(p) の凸性から、最適解は「条件なしでの最適解」と「g(x) = K の場合の最適解」のいずれか
*/
//
// verified
//   ABC 218 H - Red and Blue Lamps
//     https://atcoder.jp/contests/abc218/tasks/abc218_h
//


#include <bits/stdc++.h>
using namespace std;


/*
    min_{x}: f(x) s.t. g(x) = K
    = max_{λ は整数}: min_{x} (f(x) + λ(g(x) - K))

    max_{x}: f(x) s.t. g(x) = K
    = min_{λ は整数}: max_{x} (f(x) + λ(g(x) - K))

    【仮定】
　　・f(x), g(x): 0 以上の整数値をとる関数
　　・K: 整数
　　・g(x) = p となる x についての f(x) の最小値を h(p) とおくと、h(p) は p について下に凸な関数
　　・h(K) < ∞ (g(x) = K となる x が存在)

    【インターフェース】
　　・lag(λ): x✳︎ = argmin_{x} {f(x) + λ(g(x) - K)} としたときのペア値 {f(x✳︎), g(x✳︎)} を返す関数
　　・返り値: {最適値, そのときの λ}
*/
template<class VAL, class FUNC> pair<VAL, VAL> AliensTrick(const FUNC &lag, VAL K) {
    // zero check
    auto [f0, g0] = lag(VAL(0));
    if (g0 == K) return {f0, VAL(0)};

    // find an initial value of binary search (by doubling)
    VAL lo = 0, hi = 0, flo = f0, fhi = f0, glo = g0, ghi = g0;
    if (g0 < K) {
        lo = VAL(-1);
        tie(flo, glo) = lag(lo);
        while (glo < K) {
            assert(lo < hi);
            VAL dif = hi - lo;
            hi = lo, fhi = flo, ghi = glo, lo -= VAL(2) * dif;
            tie(flo, glo) = lag(lo);
        }
    } else {
        hi = VAL(1);
        tie(fhi, ghi) = lag(hi);
        while (ghi > K) {
            assert(lo < hi);
            VAL dif = hi - lo;
            lo = hi, flo = fhi, glo = ghi, hi += VAL(2) * dif;
            tie(fhi, ghi) = lag(hi);
        }
    }
    if (glo == K) return {flo, lo};
    if (ghi == K) return {fhi, hi};
    assert(glo > K && ghi < K);

    // binary search
    while (hi - lo > VAL(1)) {
        VAL mid = (lo + hi) / 2;
        auto [fmid, gmid] = lag(mid);
        if (gmid == K) return {fmid, mid};
        if (gmid > K) lo = mid, flo = fmid, glo = gmid;
        else hi = mid, fhi = fmid, ghi = gmid;
    }

    // in case: h(p) has a collinear part (= could not achieve g(x) = K)
    // return the Lagrange value
    VAL reslo = flo + lo * (glo - K), reshi = fhi + hi * (ghi - K);
    if (reslo >= reshi) return {reslo, lo};
    else return {reshi, hi};
}
template<class VAL, class FUNC> pair<VAL, VAL> AliensTrickMinimize(const FUNC &lag, VAL K) {
    return AliensTrick(lag, K);
}
template<class VAL, class FUNC> pair<VAL, VAL> AliensTrickMaximize(const FUNC &lag, VAL K) {
    auto reverse_lag = [&](VAL lambda) -> pair<VAL, VAL> {
        auto [f, g] = lag(-lambda);
        return {-f, g};
    };
    auto [f, g] = AliensTrick(reverse_lag, K);
    return {-f, -g};
}


//------------------------------//
// Examples
//------------------------------//

// ABC 218 H - Red and Blue Lamps
template<class VAL> struct MongeShortestPath {
    VAL INF = numeric_limits<VAL>::max() / 2;
    int CNT_INF = numeric_limits<int>::max() / 2;

    // results
    vector<VAL> dp;
    vector<int> cnt, prev;

    // solver
    template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f, bool minimize_cnt = true) {
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;

        auto relax = [&](int l, int r) -> void {
            VAL val = dp[l] + f(l, r);
            int c = cnt[l] + 1;
            if (dp[r] > val || (dp[r] == val && (minimize_cnt ? c < cnt[r] : c > cnt[r]))) {
                dp[r] = val;
                cnt[r] = c;
                prev[r] = l;
            }
        };
        auto rec = [&](auto &&rec, int l, int r) -> void {
            if (r - l <= 1) return;
            int m = (l + r) / 2;
            for (int k = prev[l]; k <= prev[r]; k++) relax(k, m);
            rec(rec, l, m);
            for (int k = l + 1; k <= m; k++) relax(k, r);
            rec(rec, m, r);
        };

        if (N > 0) relax(0, N), rec(rec, 0, N);
        vector<pair<VAL, int>> res(N + 1, make_pair(numeric_limits<VAL>::max() / 2, -1));
        res[0].first = VAL(0);
        for (int i = 1; i <= N; i++) res[i] = {dp[i], prev[i]};
        return res;
    }

    vector<int> reconstruct() {
        int N = (int)dp.size() - 1;
        vector<int> path;
        for (int v = N; v > 0; v = prev[v]) path.emplace_back(v);
        path.emplace_back(0);
        reverse(path.begin(), path.end());
        return path;
    }
};
template<class VAL> struct AliensDP {
    // inner values
    MongeShortestPath<VAL> msp;

    // solver
    template<class FUNC> pair<VAL, VAL> solve(int N, const FUNC &f, VAL D) {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            const auto &obj = msp.solve(N, fg);
            auto cnt = msp.cnt[N], res = obj.back().first - lambda * cnt;
            return {res, cnt};
        };
        return AliensTrick(lag, D);
    };

    // debugger
    template<class FUNC> void debug(int N, const FUNC &f, VAL start = 0, VAL goal = 10) {
        auto lag = [&](VAL lambda) -> void {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            const auto &dp = msp.solve(N, fg);
            auto cnt = msp.cnt[N], res = dp.back().first - lambda * cnt;
            cout << lambda << ": " << res << "(" << cnt << "); (";
            for (int i = 0; i < (int)dp.size(); i++) {
                if (i) cout << " ";
                cout << "<" << dp[i].first << ", " << dp[i].second << ">";
            }
            cout << ")" << endl;
        };
        for (VAL lambda = start; lambda <= goal; lambda++) lag(lambda);
    }
};
void ABC_218_H() {
    long long N, R;
    cin >> N >> R;
    vector<long long> A(N - 1);
    for (int i = 0; i < N - 1; i++) cin >> A[i];
    auto cost = [&](int i, int j) -> long long {
        if (j - i <= 1) return 0;
        if (i == 0) return -A[j - 2];
        else return -A[i - 1] - A[j - 2];
    };
    AliensDP<long long> ad;
    //ad.debug(N, cost);
    auto [resR, cntR] = ad.solve(N, cost, R);
    auto [resB, cntB] = ad.solve(N, cost, N - R);
    auto res = max(-resR, -resB);
    cout << res << endl;
}


int main() {
    ABC_218_H();
}