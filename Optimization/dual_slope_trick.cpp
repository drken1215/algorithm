//
// Dual Slope Trick
//
// reference:
//   maspy: slope trick (3) slope trick の凸共役
//     https://maspypy.com/slope-trick-3-slope-trick-%E3%81%AE%E5%87%B8%E5%85%B1%E5%BD%B9
//
// verified
//   yukicoder No.2114 01 Matching
//     https://yukicoder.me/problems/no/2114
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


// Dual Slope Trick
/*
    f(x): 区分線形凸関数
    ・f0: f(0)
    ・offset_l, offset_r: L, R 全体に対して加算する値
    ・L, R: x <= 0, x >= 0 の領域での線分の傾きの列
*/
template<class COORD, class SLOPE> struct DualSlopeTrick {
    using POINT = pair<COORD, SLOPE>;
    using LINE = vector<pair<COORD, SLOPE>>;

    // inner data
    SLOPE f0, offsetL, offsetR, INF;
    priority_queue<SLOPE> L;
    priority_queue<SLOPE, vector<SLOPE>, greater<SLOPE>> R;

    // constructors
    // initialize: f(x) = 0 (x == 0), INF (otherwise)
    DualSlopeTrick(SLOPE inf = numeric_limits<SLOPE>::max() / 2)
        : f0(0), offsetL(0), offsetR(0), INF(inf) {
        assert(inf > 0);
    }
    DualSlopeTrick(const DualSlopeTrick&) = default;
    DualSlopeTrick& operator = (const DualSlopeTrick&) = default;

    // basic operations
    constexpr int sizeL() const { return (int)L.size(); }
    constexpr int sizeR() const { return (int)R.size(); }
    constexpr int size() const { return sizeL() + sizeR(); }
    constexpr void pushL(const SLOPE &v) { L.push(v - offsetL); }
    constexpr void pushR(const SLOPE &v) { R.push(v - offsetR); }
    constexpr SLOPE topL() const { return L.empty() ? -INF : L.top() + offsetL; }
    constexpr SLOPE topR() const { return R.empty() ? INF : R.top() + offsetR; }
    constexpr SLOPE popL() {
        auto res = topL();
        if (!L.empty()) L.pop();
        return res;
    }
    constexpr SLOPE popR() {
        auto res = topR();
        if (!R.empty()) R.pop();
        return res;
    }

    // getter and debugger
    constexpr SLOPE get_f0() const { return f0; }  // O(1)
    constexpr pair<LINE, LINE> get_lines() const {  // O(N log N)
        auto L2 = L;
        auto R2 = R;
        COORD lx = 0, ly = f0, rx = 0, ry = f0;
        LINE resL{{lx, ly}}, resR{{rx, ry}};
        while (!L2.empty()) {
            auto dif = L2.top() + offsetL;
            L2.pop();
            lx--, ly -= dif;
            resL.emplace_back(lx, ly);
        }
        while (!R2.empty()) {
            auto dif = R2.top() + offsetR;
            R2.pop();
            rx++, ry += dif;
            resR.emplace_back(rx, ry);
        }
        return {resL, resR};
    }
    constexpr SLOPE eval(COORD x) const {  // O(N log N)
        SLOPE res = f0;
        auto L2 = L;
        auto R2 = R;
        if (x > L2.size() || x < -R2.size()) return INF;
        if (x >= 0) {
            for (int i = 0; i < x; i++) {
                auto dif = L2.top() + offsetL;
                L2.pop();
                res -= dif;
            }
        } else {
            for (int i = 0; i < -x; i++) {
                auto dif = R2.top() + offsetR;
                R2.pop();
                res += dif;
            }
        }
        return res;
    }
    constexpr SLOPE get_min() const {  // O(N log N)
        SLOPE res = f0;
        auto L2 = L;
        auto R2 = R;
        while (!L2.empty()) {
            auto dif = L2.top() + offsetL;
            L2.pop();
            res += min(SLOPE(0), -dif);
        }
        while (!R2.empty()) {
            auto dif = R2.top() + offsetR;
            R2.pop();
            res += min(SLOPE(0), dif);
        }
        return res;
    }
    constexpr friend ostream &operator << (ostream &os, DualSlopeTrick st) {  // O(N log N)
        auto [lineL, lineR] = st.get_lines();
        os << endl << "left: ";
        for (auto [x, y] : lineL) os << "(" << x << ", " << y << ") ";
        os << endl << "right: ";
        for (auto [x, y] : lineR) os << "(" << x << ", " << y << ") ";
        return os << endl;
    }
    
    // f(x) += b, O(1)
    DualSlopeTrick &add_const(const SLOPE &b) {
        f0 += b;
        return *this;
    }

    // f(x) += ax + b, O(1)
    DualSlopeTrick &add_linear(const SLOPE &a, const SLOPE &b) {
        offsetL += a, offsetR += a;
        return add_const(b);
    }

    // f(x) += max(0, c(x - a)), O(|a| log N)
    DualSlopeTrick &add_relu(const SLOPE &c, COORD a) {
        slide(-a);
        if (c > SLOPE(0)) offsetR += c;
        else offsetL += c;
        slide(a);
        return *this;
    }
    DualSlopeTrick &add_relu(COORD a) {
        return add_relu(1, a);
    }
    DualSlopeTrick &add_irelu(COORD a) {
        return add_relu(-1, a);
    }

    // f(x) += c|x - a|, O(|a| log N)
    DualSlopeTrick &add_abs(const SLOPE &c, COORD a) {
        add_relu(c, a), add_relu(-c, a);
        return *this;
    }
    DualSlopeTrick &add_abs(COORD a) {
        return add_abs(1, a);
    }

    // f(x) <- g(x) = f(x - 1), O(log N)
    DualSlopeTrick &slide() {
        SLOPE X = popL();
        pushR(X);
        return add_const(-X);
    }

    // f(x) <- g(x) = f(x + 1), O(log N)
    DualSlopeTrick &slide_rev() {
        SLOPE X = popR();
        pushL(X);
        return add_const(X);
    }

    // f(x) <- g(x) = f(x - a), O(|a| log N)
    DualSlopeTrick &slide(COORD a) {
        while (a > 0) a--, slide();
        while (a < 0) a++, slide_rev();
        return *this;
    }

    // f(x) <- g(x) = min_{a <= y <= b} f(x - y) = min_{x-b <= y <= x-a} f(y), O((|a| + |b|) log N)
    DualSlopeTrick &slide(COORD a, COORD b) {
        assert(a <= b);
        slide(a);
        for (COORD i = 0; i < b - a; i++) {
            add_const(min(SLOPE(0), -topL()));
            pushL(0), pushR(popL());
        }
        return *this;
    }
};


//------------------------------//
// Examples
//------------------------------//

// yukicoder No.2114 01 Matching
/*
    v1, ..., vN: 座標を小さい順に並べたもの。赤か青。赤を source、青を sink にする
    dp_{i}[x] := 最初の i 個を終えた時点で染み出しフローが x であるときの、染み出し分も含めたコストの最小値

    ・i 番目が赤のとき（D = x[i+1] - x[i] とする）
    nex[x] = dp[x-1] + D|x|（平行移動 + 直線加算）

    ・i 番目が青のとき（D = x[i+1] - x[i] とする）
    nex[x] = min(dp[x], dp[x+1]) + D|x|（スライド最小値 + 直線加算）
*/
void yukicoder_2114() {
    using i128 = __int128_t;
    const int RED = 0, BLUE = 1;
    long long N, M, K;
    cin >> N >> M >> K;
    vector<long long> B(N), R(M);
    map<long long, vector<pair<long long, int>>> mp;
    for (int i = 0; i < N; i++) cin >> B[i];
    for (int i = 0; i < M; i++) cin >> R[i];
    if (N < M) swap(N, M), swap(B, R);  // RED の方が少なくする
    for (int i = 0; i < N; i++) mp[B[i] % K].emplace_back(B[i] / K, BLUE);
    for (int i = 0; i < M; i++) mp[R[i] % K].emplace_back(R[i] / K, RED);
    
    auto calc = [&](vector<pair<long long, int>> &v) -> pair<bool, long long> {
        sort(v.begin(), v.end());
        DualSlopeTrick<long long, i128> st(1LL<<50);  // INF を積み上げても overflow しないように！
        for (int i = 0; i < v.size(); i++) {
            long long D = (i + 1 < v.size() ? v[i + 1].first - v[i].first : 0);
            if (v[i].second == RED) {
                st.slide(1);
                st.add_abs(D, 0);
            } else {
                st.slide(-1, 0);
                st.add_abs(D, 0);
            }
        }
        auto res = st.get_f0();
        return {bool(res < st.INF/2), res};
    };
    bool can = true;
    long long res = 0;
    for (auto [key, v] : mp) {
        auto [feasible, tmp] = calc(v);
        if (!feasible) can = false;
        else res += tmp;
    }
    cout << (can ? res : -1) << endl;
}


int main() {
    yukicoder_2114();
}