//
// Slope Trick
//
// verified
//   第2回 ドワンゴからの挑戦状 予選 E - 花火
//     https://atcoder.jp/contests/dwango2016-prelims/tasks/dwango2016qual_e
//
//   AtCoder ABC 217 H - Snuketoon
//     https://atcoder.jp/contests/abc217/tasks/abc217_h
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
    constexpr pair<LINE, LINE> get_lines() {
        LINE resL, resR;
        vector<COORD> L2, R2;
        COORD sumL = 0, sumR = 0;
        int ln = 0, rn = 0;
        while (!L.empty()) {
            auto x = popL();
            sumL += x;
            auto y = min_f + sumL - x * (++ln);
            if (resL.empty() || x != resL.back().first) resL.emplace_back(x, y);
            L2.emplace_back(x);
        }
        while (!R.empty()) {
            auto x = popR();
            sumR += x;
            auto y = min_f + x * (++rn) - sumR;
            if (resR.empty() || x != resR.back().first) resR.emplace_back(x, y);
            R2.emplace_back(x);
        }
        for (auto v : L2) L.push(v);
        for (auto v : R2) R.push(v);
        return {resL, resR};
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

    // f(x) <- g(x) = f(x - a), O(1)
    SlopeTrick &slide(const COORD &a) {
        offsetL += a, offsetR += a;
        return *this;
    }

    // f(x) <- g(x) = min_{a <= y <= b} f(x - y) = min_{x-b <= y <= x-a} f(y), O(1)
    SlopeTrick &slide(const COORD &a, const COORD &b) {
        assert(a <= b);
        offsetL += a, offsetR += b;
        return *this;
    }

    // f(x) <- g(x) = min_{0 <= y <= a} f(x - y) = min_{x-a <= y <= x} f(y), O(1)
    SlopeTrick &slide_right_curve_to_right(const COORD &a) {
        assert(a >= 0);
        offsetR += a;
        return *this;
    }

    // f(x) <- g(x) = min_{0 <= y <= a} f(x + y) = min_{x <= y <= x+a} f(y),  O(1)
    SlopeTrick &slide_left_curve_to_left(const COORD &a) {
        assert(a >= 0);
        offsetL -= a;
        return *this;
    }
};


//------------------------------//
// Solver
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


int main() {
    //DWANGO_2nd_prelims_E();
    ABC_217_H();
}