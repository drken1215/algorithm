//
// Slope Trick
//
// reference:
//   maspy: slope trick (1) 解説編
//     https://maspypy.com/slope-trick-1-%E8%A7%A3%E8%AA%AC%E7%B7%A8
//
// verified
//   第2回 ドワンゴからの挑戦状 予選 E - 花火 (for 累積min, +|x-a|)
//     https://atcoder.jp/contests/dwango2016-prelims/tasks/dwango2016qual_e
//
//   AtCoder AWC 0100 N - 株価の補正 (for slide(1, INF) など)
//     https://atcoder.jp/contests/awc0100/tasks/awc0100_n
//
//   AtCoder ABC 217 H - Snuketoon (for slide(-dt, dt), +max(0, a-x) など)
//     https://atcoder.jp/contests/abc217/tasks/abc217_h
//
//   AtCoder ARC 123 D - Inc, Dec - Decomposition
//     https://atcoder.jp/contests/arc123/tasks/arc123_d
//
//   KUPC 2016 H - 壁壁壁壁壁壁壁 (for slide(-INF, A[i]-B[i]), eval(x) など)
//     https://atcoder.jp/contests/kupc2016/tasks/kupc2016_h
//
//   UTPC 2012 L - じょうしょうツリー (for add(g))
//     https://atcoder.jp/contests/utpc2012/tasks/utpc2012_12
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


// Slope Trick
/*
    f(x): 区分線形凸関数
    ・min_f: f(x) の最小値
    ・offset_l, offset_r: L, R 全体に対して加算する値
    ・L, R: f(x) の傾きが 1 変化する x 座標の多重集合 (min_f の左側と右側)
*/
template<class COORD> struct SlopeTrick {
    using POINT = pair<COORD, COORD>;
    using LINE = vector<pair<COORD, COORD>>;

    // inner data
    COORD min_f, offsetL, offsetR, INF;
    priority_queue<COORD> L;
    priority_queue<COORD, vector<COORD>, greater<COORD>> R;

    // constructors
    SlopeTrick(COORD inf = numeric_limits<COORD>::max() / 2)
        : min_f(0), offsetL(0), offsetR(0), INF(inf) {
        assert(inf > 0);
    }
    SlopeTrick(const SlopeTrick&) = default;
    SlopeTrick& operator = (const SlopeTrick&) = default;

    // getter and debugger
    constexpr int sizeL() const { return (int)L.size(); }
    constexpr int sizeR() const { return (int)R.size(); }
    constexpr int size() const { return sizeL() + sizeR(); }
    constexpr void pushL(const COORD &v) { L.push(v - offsetL); }
    constexpr void pushR(const COORD &v) { R.push(v - offsetR); }
    constexpr COORD topL() const { return L.empty() ? -INF : L.top() + offsetL; }
    constexpr COORD topR() const { return R.empty() ? INF : R.top() + offsetR; }
    constexpr COORD popL() {
        auto res = topL();
        if (!L.empty()) L.pop();
        return res;
    }
    constexpr COORD popR() {
        auto res = topR();
        if (!R.empty()) R.pop();
        return res;
    }
    constexpr COORD get_min() const { return min_f; }
    constexpr pair<COORD, COORD> get_argmin() const { return {topL(), topR()}; }
    constexpr pair<LINE, LINE> get_lines() const {
        LINE resL, resR;
        auto L2 = L;
        auto R2 = R;
        COORD sumL = 0, sumR = 0;
        int ln = 0, rn = 0;
        while (!L2.empty()) {
            auto x = L2.top() + offsetL;
            L2.pop();
            sumL += x;
            auto y = min_f + sumL - x * (++ln);
            if (resL.empty() || x != resL.back().first) resL.emplace_back(x, y);
        }
        while (!R2.empty()) {
            auto x = R2.top() + offsetR;
            R2.pop();
            sumR += x;
            auto y = min_f + x * (++rn) - sumR;
            if (resR.empty() || x != resR.back().first) resR.emplace_back(x, y);
        }
        return {resL, resR};
    }
    constexpr COORD eval(const COORD &x) const {
        COORD res = 0;
        auto L2 = L;
        auto R2 = R;
        while (!L2.empty()) {
            auto t = L2.top() + offsetL;
            L2.pop();
            res += max(COORD(0), t - x);
        }
        while (!R2.empty()) {
            auto t = R2.top() + offsetR;
            R2.pop();
            res += max(COORD(0), x - t);
        }
        return res + min_f;
    }
    constexpr friend ostream &operator << (ostream &os, SlopeTrick st) {
        auto [lineL, lineR] = st.get_lines();
        os << endl << "left: ";
        for (auto [x, y] : lineL) os << "(" << x << ", " << y << ") ";
        os << endl << "right: ";
        for (auto [x, y] : lineR) os << "(" << x << ", " << y << ") ";
        return os << endl;
    }
    
    // f(x) += b, O(1)
    SlopeTrick &add_const(const COORD &b) {
        min_f += b;
        return *this;
    }

    // f(x) += max(0, x - a), O(log N)
    SlopeTrick &add_relu(const COORD &a) {
        if (topL() <= a) pushR(a);
        else {
            min_f += max(COORD(0), topL() - a);
            pushL(a), pushR(popL());
        }
        return *this;
    }

    // f(x) += max(0, a - x), O(log N)
    SlopeTrick &add_irelu(const COORD &a) {
        if (topR() >= a) pushL(a);
        else {
            min_f += max(COORD(0), a - topR());
            pushR(a), pushL(popR());
        }
        return *this;
    }

    // f(x) += |x - a|, O(log N)
    SlopeTrick &add_abs(const COORD &a) {
        add_relu(a), add_irelu(a);
        return *this;
    }

    // f(x) <- g(x) = min_{y <= x} f(y), (\_/ -> \__), O(1)
    SlopeTrick &clear_right() {
        R = priority_queue<COORD, vector<COORD>, greater<COORD>>{};
        return *this;
    }

    // f(x) <- g(x) = min_{y >= x} f(y), (\_/ -> __/), O(1)
    SlopeTrick &clear_left() {
        L = priority_queue<COORD>{};
        return *this;
    }

    // f(x) <- 0
    SlopeTrick &clear() {
        *this = SlopeTrick();
        return *this;
    }

    // f(x) <- g(x) = f(x - a), O(1)
    SlopeTrick &slide(const COORD &a) {
        offsetL += a, offsetR += a;
        return *this;
    }

    // f(x) <- g(x) = min_{a <= y <= b} f(x - y) = min_{x-b <= y <= x-a} f(y), O(1)
    SlopeTrick &slide(const COORD &a, const COORD &b) {
        assert(a <= b);
        if (a <= -INF) clear_left();
        else offsetL += a;
        if (b >= INF) clear_right();
        else offsetR += b;
        return *this;
    }

    // f(x) <- g(x) = min_{0 <= y <= a} f(x - y) = min_{x-a <= y <= x} f(y), O(1)
    SlopeTrick &slide_right_curve_to_right(const COORD &a) {
        assert(a >= 0);
        slide(0, a);
        return *this;
    }

    // f(x) <- g(x) = min_{-a <= y <= 0} f(x - y) = min_{x <= y <= x+a} f(y),  O(1)
    SlopeTrick &slide_left_curve_to_left(const COORD &a) {
        assert(a >= 0);
        slide(-a, 0);
        return *this;
    }

    // f(x) += g(x), O((log N)^2)
    // attention: this function is destructive
    SlopeTrick &add(SlopeTrick &g) {
        if (size() < g.size()) {
            swap(min_f, g.min_f);
            swap(offsetL, g.offsetL), swap(offsetR, g.offsetR);
            swap(L, g.L), swap(R, g.R);
        }
        min_f += g.min_f;
        while (g.L.size()) add_irelu(g.popL());
        while (g.R.size()) add_relu(g.popR());
        return *this;
    }
};


//------------------------------//
// Examples
//------------------------------//

// 第2回 ドワンゴからの挑戦状 予選 E - 花火
/*
    dp[x] := 今考えている花火までについての最終位置が x のときのスコア
    nex[x] = min_{y <= x} (dp[y]) + |x - v[0]| + |x - v[1]] + ... + |x - v[K-1]|

    min_{y <= x} (dp[y]) は clear_right()
*/
void DWANGO_2nd_prelims_E() {
    long long N, L, INF = 1LL<<50;
    cin >> N >> L;
    map<long long, vector<long long>> mp;
    for (int i = 0; i < N; i++) {
        long long t, P;
        cin >> t >> P;
        mp[t].emplace_back(P);
    }
    SlopeTrick<long long> dp;
    for (auto [key, v] : mp) {
        dp.clear_right();
        for (auto a : v) dp.add_abs(a);
    }
    long long res = dp.min_f;
    cout << res << endl;
}


// AtCoder AWC 0100 N - 株価の補正
/*
    nex[x] = min_{y <= x-1} dp[y] + |x - H[i]|
*/
void AWC_0100_N() {
    long long N;
    cin >> N;
    vector<long long> H(N);
    for (int i = 0; i < N; i++) cin >> H[i];

    SlopeTrick<long long> st;
    for (int i = 0; i < N; i++) {
        st.slide(1, st.INF);
        st.add_abs(H[i]);
    }
    cout << st.get_min() << endl;
}


// AtCoder ABC 217 H - Snuketoon
/*
    dp[x] := 最後の攻撃の時に位置 x にいる場合の、これまでの総ダメージの最小値

    Di = 0 のとき
        nex[x] = min_{x-dt <= y <= x+dt} dp[y] + max(0, Xi - x)
    Di = 1 のとき
        nex[x] = min_{x-dt <= y <= x+dt} dp[y] + max(0, x - Xi)
*/
void ABC_217_H() {
    int N;
    cin >> N;
    vector<long long> T(N), D(N), X(N);
    for (int i = 0; i < N; i++) cin >> T[i] >> D[i] >> X[i];
    SlopeTrick<long long> dp;
    for (int i = 0; i < N * 5; i++) dp.add_abs(0);  // 擬似的に x = 0 以外を ∞ にする
    long long prev = 0;
    for (int i = 0; i < N; i++) {
        long long dt = T[i] - prev;
        prev = T[i];
        dp.slide(-dt, dt);
        if (D[i] == 0) dp.add_irelu(X[i]);
        else dp.add_relu(X[i]);
    }
    long long res = dp.min_f;
    cout << res << endl;
}


// AtCoder ARC 123 D - Inc, Dec - Decomposition
/*
    dp[x] := 最初の i 項目までについて、正数列 B の末項が x である場合の操作回数の最小値
        y = B[i-1], x = B[i] を満たす条件は
        ・y <= x
        ・A[i-1] - y >= A[i] - x 　⇔  y <= x - (A[i] - A[i-1])
    nex[x] = min_{y <= x - max(0, A[i]-A[i-1])} dp[y] + |x| + |x - A[i]|
*/
void ARC_123_D() {
    long long N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    SlopeTrick<long long> dp;
    dp.add_abs(0), dp.add_abs(A[0]);
    for (int i = 1; i < N; i++) {
        dp.slide(max(0LL, A[i]-A[i-1]), dp.INF);
        dp.add_abs(0), dp.add_abs(A[i]);
    }
    auto res = dp.get_min();
    cout << res << endl;
}


// KUPC 2016 H - 壁壁壁壁壁壁壁
/*
    dp[x] := その時点で、はみ出し枚数が x のときのコストの最小値 (最後のはみ出しを含む)
    nex[x] = min_{x - (A[i] - B[i]) <= y} (dp[y]) + |x|
*/
void KUPC_2016_H() {
    long long N;
    cin >> N;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> B[i];
    SlopeTrick<long long> dp;
    for (int i = 0; i < N * 2; i++) dp.add_abs(0);
    for (int i = 0; i < N; i++) {
        dp.slide(-dp.INF, A[i]-B[i]);
        dp.add_abs(0);
    }
    long long res = dp.eval(0);
    cout << res << endl;
}


// UTPC 2012 L - じょうしょうツリー
/*
    深さを足しておくことで、条件を x[p] >= x[c] と表せる状態にしておく

    dp_v[x] = Σ_{c}(min_{y <= x}(dp_c[y])) + |x - A[v]|
*/
void UTPC_2012_L() {
    int N;
    cin >> N;
    vector<vector<int>> tree(N);
    vector<long long> P(N, -1), A(N);
    cin >> A[0];
    for (int i = 1; i < N; i++) {
        cin >> P[i] >> A[i], P[i]--;
        tree[P[i]].emplace_back(i);
    }
    vector<SlopeTrick<long long>> dp(N);
    auto rec = [&](auto &&rec, int v, int depth = 0) -> void {
        A[v] += depth;
        for (auto c : tree[v]) {
            rec(rec, c, depth + 1);
            dp[c].clear_right();
            dp[v].add(dp[c]);
        }
        dp[v].add_abs(A[v]);
    };
    rec(rec, 0);
    cout << dp[0].get_min() << endl;
}


int main() {
    //DWANGO_2nd_prelims_E();
    //AWC_0100_N();
    //ABC_217_H();
    //ARC_123_D();
    //KUPC_2016_H();
    UTPC_2012_L();
}