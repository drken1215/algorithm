//
// Monge 単一始点 d-辺最短路問題 (by Alien DP), O(N (log N)^2)
//   最短路問題を解くのには noshi's 簡易 LARSCH を用いる
//
// verified
//   ABC 218 H - Red and Blue Lamps
//     https://atcoder.jp/contests/abc218/tasks/abc218_h
//


#include <bits/stdc++.h>
using namespace std;


/*
    Alien's Trick (by ラグランジュ緩和)

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
template<class VAL> struct AliensTrick {
    // find lo such that g(lo) >= R
    template<class FUNC> tuple<VAL, VAL, VAL> find_low(const FUNC &lag, VAL R) {
        auto [f0, g0] = lag(VAL(0));
        VAL lo = 0, hi = 0, flo = f0, glo = g0;
        if (glo >= R) return tie(lo, flo, glo);
        lo = VAL(-1);
        tie(flo, glo) = lag(lo);
        while (glo < R) {
            assert(lo < hi);
            VAL dif = hi - lo;
            hi = lo, lo -= VAL(2) * dif;
            tie(flo, glo) = lag(lo);
        }
        return tie(lo, flo, glo);
    }

    // find hi such that g(hi) <= L
    template<class FUNC> tuple<VAL, VAL, VAL> find_high(const FUNC &lag, VAL L) {
        auto [f0, g0] = lag(VAL(0));
        VAL lo = 0, hi = 0, fhi = f0, ghi = g0;
        if (ghi <= L) return tie(hi, fhi, ghi);
        hi = VAL(1);
        tie(fhi, ghi) = lag(hi);
        while (ghi > L) {
            assert(lo < hi);
            VAL dif = hi - lo;
            lo = hi, hi += VAL(2) * dif;
            tie(fhi, ghi) = lag(hi);
        }
        return tie(hi, fhi, ghi);
    }

    // solver
    template<class FUNC> pair<VAL, VAL> solve(const FUNC &lag, VAL K) {
        // zero check
        auto [f0, g0] = lag(VAL(0));
        if (g0 == K) return {f0, VAL(0)};

        // find an initial value of binary search (by doubling)
        auto [lo, flo, glo] = find_low(lag, K);
        auto [hi, fhi, ghi] = find_high(lag, K);
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
    template<class FUNC> pair<VAL, VAL> minimize(const FUNC &lag, VAL K) {
        return solve(lag, K);
    }
    template<class FUNC> pair<VAL, VAL> maximize(const FUNC &lag, VAL K) {
        auto reverse_lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto [f, g] = lag(-lambda);
            return {-f, g};
        };
        auto [f, g] = solve(reverse_lag, K);
        return {-f, -g};
    }

    // in case: L <= g(x) <= R
    template<class FUNC> pair<VAL, VAL> solve(const FUNC &lag, VAL L, VAL R) {
        assert(L <= R);
        auto [f0, g0] = lag(VAL(0));
        if (L <= g0 && g0 <= R) return {f0, VAL(0)};
        if (g0 > R) return solve(lag, R);
        else return solve(lag, L);
    }
    template<class FUNC> pair<VAL, VAL> minimize(const FUNC &lag, VAL L, VAL R) {
        return solve(lag, L, R);
    }
    template<class FUNC> pair<VAL, VAL> maximize(const FUNC &lag, VAL L, VAL R) {
        auto reverse_lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto [f, g] = lag(-lambda);
            return {-f, g};
        };
        auto [f, g] = solve(reverse_lag, L, R);
        return {-f, -g};
    }
};

// noshi's simplified LARSCH
// find the shortest path from vertex 0 on DAG in O(N log N)
// vertex: 0, 1, 2, ..., N, f(i, j) must be Monge
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

// Alien DP
// find the d-edges shortest path from vertex 0 on DAG by Lagrange relaxation
// vertex: 0, 1, 2, ..., N, f(i, j) must be Monge
template<class VAL> struct AlienDP {
    // inner values
    MongeShortestPath<VAL> msp;
    AliensTrick<VAL> at;

    // solver
    template<class FUNC> pair<VAL, VAL> solve(int N, const FUNC &f, VAL D) {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            const auto &obj = msp.solve(N, fg);
            auto cnt = msp.cnt[N];
            auto res = obj.back().first - lambda * cnt;
            return {res, cnt};
        };
        return at.solve(lag, D);
    };

    // debugger
    template<class FUNC> void debug(int N, const FUNC &f, VAL start = 0, VAL goal = 10) {
        auto lag = [&](VAL lambda) -> void {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            const auto &dp = msp.solve(N, fg);
            auto cnt = msp.cnt[N];
            auto res = dp.back().first - lambda * cnt;
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


//------------------------------//
// Examples
//------------------------------//

// ABC 218 H - Red and Blue Lamps
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
    AlienDP<long long> ad;
    //ad.debug(N, cost);
    auto [resR, cntR] = ad.solve(N, cost, R);
    auto [resB, cntB] = ad.solve(N, cost, N - R);
    auto res = max(-resR, -resB);
    cout << res << endl;
}


int main() {
    ABC_218_H();
}