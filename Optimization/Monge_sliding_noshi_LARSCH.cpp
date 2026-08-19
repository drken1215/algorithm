//
// noshi's 簡易 LARSCH スライドアクセス対応版
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


// noshi's simplified LARSCH
// find shortest path from vertex 0 on DAG in O(N log N)
// vertex: 0, 1, 2, ..., N, f(i, j) must be Monge
template<class VAL> struct MongeShortestPath {
    VAL INF = numeric_limits<VAL>::max() / 2;
    int CNT_INF = numeric_limits<int>::max() / 2;

    // results
    vector<VAL> dp;
    vector<int> cnt, prev;

    // solver (random access ver)
    template<class FUNC> vector<pair<VAL, int>> solve(int N, const FUNC &f, bool minimize_cnt = true) {
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;

        auto relax = [&](int l, int r) -> void {
            VAL val = dp[l] + f(l, r);
            int cn = cnt[l] + 1;
            if (dp[r] > val || (dp[r] == val && (minimize_cnt ? cn < cnt[r] : cn > cnt[r]))) {
                dp[r] = val;
                cnt[r] = cn;
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

    // solver (slide access ver)
    template<class STATE, class ADD, class DEL, class GETCOST> vector<pair<VAL, int>> solve
    (int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, bool minimize_cnt = true) {
        dp.assign(N + 1, INF);
        cnt.assign(N + 1, CNT_INF);
        prev.assign(N + 1, 0);
        dp[0] = 0, cnt[0] = 0;

        STATE win[2] = { ini, ini };
        int cl[2] = {0, 0}, cr[2] = {0, 0};

        auto move_cursor = [&](int c, int l, int r) {
            while (cr[c] < r) add(win[c], cr[c]++);
            while (cl[c] > l) add(win[c], --cl[c]);
            while (cr[c] > r) del(win[c], --cr[c]);
            while (cl[c] < l) del(win[c], cl[c]++);
        };
        auto relax = [&](int c, int l, int r) {
            move_cursor(c, l, r);
            VAL val = dp[l] + get(win[c]);
            int cn = cnt[l] + 1;
            if (dp[r] > val || (dp[r] == val && (minimize_cnt ? cn < cnt[r] : cn > cnt[r]))) {
                dp[r] = val;
                cnt[r] = cn;
                prev[r] = l;
            }
        };
        auto rec = [&](auto &&rec, int l, int r) -> void {
            if (r - l <= 1) return;
            int m = (l + r) / 2;
            for (int k = prev[l]; k <= prev[r]; k++) relax(0, k, m);
            rec(rec, l, m);
            for (int k = l + 1; k <= m; k++) relax(1, k, r);
            rec(rec, m, r);
        };

        if (N > 0) relax(0, 0, N), rec(rec, 0, N);
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

// Alien DP
// find the d-edges shortest path from vertex 0 on DAG by Lagrange relaxation
// vertex: 0, 1, 2, ..., N, f(i, j) must be Monge
template<class VAL> struct AlienDP {
    // inner values
    MongeShortestPath<VAL> msp;
    AliensTrick<VAL> at;

    // solver
    template<class FUNC> pair<VAL, VAL> solve(int N, const FUNC &f, int D, const string &solver = "simple_larsch") {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            vector<pair<VAL, int>> obj;
            if (solver == "simple_larsch") obj = msp.solve(N, fg);
            auto cnt = msp.cnt[N];
            auto res = obj.back().first - lambda * cnt;
            return {res, cnt};
        };
        return at.solve(lag, D);
    }

    // solver (slide access ver.)
    template<class STATE, class ADD, class DEL, class GETCOST> pair<VAL, int> solve
    (int N, const STATE &ini, const ADD &add, const DEL &del, const GETCOST &get, int D) {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto get_lag = [&](STATE &st) -> VAL { return get(st) + lambda; };
            vector<pair<VAL, int>> obj = msp.solve(N, ini, add, del, get_lag);
            auto cnt = msp.cnt[N];
            auto res = obj.back().first - lambda * cnt;
            return {res, cnt};
        };
        return at.solve(lag, D);
    }

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

// Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem
// by Alien DP (not the assumed solution)
void Codeforces438_F_alien() {
    long long N, K;
    cin >> N >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];

    struct STATE {
        vector<long long> nums;
        long long cost;
        STATE(int N) : nums(N + 1, 0), cost(0) {}
    };
    STATE ini(N);
    auto add = [&](STATE &st, int i) -> void {
        st.cost += st.nums[A[i]];
        st.nums[A[i]]++;
    };
    auto del = [&](STATE &st, int i) -> void {
        st.nums[A[i]]--;
        st.cost -= st.nums[A[i]];
    };
    auto get = [&](STATE &st) -> long long {
        return st.cost;
    };

    AlienDP<long long> adp;
    auto res = adp.solve(N, ini, add, del, get, K);  // we should use slide access
    cout << res.first << endl;
}


int main() {
    Codeforces438_F_alien();
}