//
// 各辺のコストをスライドアクセスで計算する Monge 最短路
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
// Utility
//------------------------------//

using ll = long long;
using i128 = __int128_t;
using u128 = __uint128_t;
using pint = pair<int, int>;
using pll = pair<long long, long long>;
using tll = array<long long, 3>;
using fll = array<long long, 4>;
using vint = vector<int>;
using vll = vector<long long>;
using dint = deque<int>;
using dll = deque<long long>;
using vvint = vector<vector<int>>;
using vvll = vector<vector<long long>>;
using vpll = vector<pair<long long, long long>>;
template<class T> using min_priority_queue = priority_queue<T, vector<T>, greater<T>>;

template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }
template<class S, class T> inline auto maxll(S a, T b) { return max(ll(a), ll(b)); }
template<class S, class T> inline auto minll(S a, T b) { return min(ll(a), ll(b)); }
template<class T> auto max(const T &a) { return *max_element(a.begin(), a.end()); }
template<class T> auto min(const T &a) { return *min_element(a.begin(), a.end()); }
template<class T> auto argmax(const T &a) { return max_element(a.begin(), a.end()) - a.begin(); }
template<class T> auto argmin(const T &a) { return min_element(a.begin(), a.end()) - a.begin(); }
template<class T> auto accum(const vector<T> &a) { return accumulate(a.begin(), a.end(), T()); }
template<class T> auto accum(const deque<T> &a) { return accumulate(a.begin(), a.end(), T()); }

#define REP(i, a) for (long long i = 0; i < (long long)(a); i++)
#define REP2(i, a, b) for (long long i = a; i < (long long)(b); i++)
#define RREP(i, a) for (long long i = (a)-1; i >= (long long)(0); --i)
#define RREP2(i, a, b) for (long long i = (b)-1; i >= (long long)(a); --i)
#define EB emplace_back
#define PF push_front
#define PB push_back
#define MP make_pair
#define FI first
#define SE second
#define ALL(x) x.begin(), x.end()
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl

// input
template<class T> istream& operator >> (istream &is, vector<T> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, deque<T> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }
template<class T> istream& operator >> (istream &is, vector<vector<T>> &P)
{ for (int i = 0; i < (int)P.size(); ++i) cin >> P[i]; return is; }

// output
template<class S, class T> ostream& operator << (ostream &s, const pair<S, T> &P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 2> &P)
{ return s << '<' << P[0] << "," << P[1] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 3> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << '>'; }
template<class T> ostream& operator << (ostream &s, const array<T, 4> &P)
{ return s << '<' << P[0] << "," << P[1] << "," << P[2] << "," << P[3] << '>'; }
template<class T> ostream& operator << (ostream &s, const vector<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const deque<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const vector<vector<T>> &P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, const set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const multiset<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, const unordered_set<T> &P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
template<class S, class T> ostream& operator << (ostream &s, const unordered_map<S, T> &P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }
void yes(bool a) { cout << (a ? "yes" : "no") << endl; }
void YES(bool a) { cout << (a ? "YES" : "NO") << endl; }
void Yes(bool a) { cout << (a ? "Yes" : "No") << endl; }
const vector<int> DX = {1, 0, -1, 0, 1, -1, 1, -1};
const vector<int> DY = {0, 1, 0, -1, 1, -1, -1, 1};


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

// noshi's simplified LARSCH
// find shortest path from vertex 0 on DAG in O(N log N)
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
            else if (solver == "monotone") tmp = MonotoneMinima<VAL>(N + 1, N + 1, f2);
            for (int i = d; i <= N; i++) dp[i] = tmp[i].first;
            res.emplace_back(dp[N]);
        }
        return res;
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
// you can choose noshi LARSCH or LARSCH
template<class VAL> struct AlienDP {
    // inner values
    MongeShortestPath<VAL> msp;
    AliensTrick<VAL> at;

    // solver
    template<class FUNC> pair<VAL, VAL> solve(int N, const FUNC &f, VAL D, const string &solver = "noshi") {
        auto lag = [&](VAL lambda) -> pair<VAL, VAL> {
            auto fg = [&](int i, int j) -> VAL { return f(i, j) + lambda; };
            vector<pair<VAL, int>> obj;
            if (solver == "noshi") obj = msp.solve(N, fg);
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
    auto res = msp.solve(N, f, K, "monotone");
    cout << res[K] << endl;
}

// Codeforces Round 438 (Div. 1 + Div. 2 combined) F. Yet Another Minimization Problem
// by Alien DP
void Codeforces438_F_alien() {
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
    AlienDP<long long> msp;
    auto res = msp.solve(N, f, K);
    cout << res.first << endl;
}


int main() {
    //Codeforces438_F_enum();
    Codeforces438_F_alien();
}