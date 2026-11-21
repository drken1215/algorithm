#pragma once
#include <bits/stdc++.h>
#include <cassert>
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
        if (x < -L2.size() || x > R2.size()) return INF;
        if (x < 0) {
            for (int i = 0; i < -x; i++) {
                auto dif = L2.top() + offsetL;
                L2.pop();
                res -= dif;
            }
        } else {
            for (int i = 0; i < x; i++) {
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

    // restrict f(x)-domain into x >= 0, O(1)
    DualSlopeTrick &clear_left() {
        L = priority_queue<COORD>{};
        offsetL = 0;
        return *this;
    }

    // restrict f(x)-domain into x <= 0, O(1)
    DualSlopeTrick &clear_right() {
        R = priority_queue<COORD, vector<COORD>, greater<COORD>>{};
        offsetR = 0;
        return *this;
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

    // f(x) <- h(x) = min_{x = y + z, z >= 0} (f(y) + g(z)), amotized O((log N)^2)
    DualSlopeTrick &min_plus_convolution(DualSlopeTrick &g) {
        if (size() < g.size()) {
            swap(f0, g.f0);
            swap(offsetL, g.offsetL), swap(offsetR, g.offsetR);
            swap(L, g.L), swap(R, g.R);
        }
        f0 += g.f0;
        while (g.L.size()) pushL(g.popL());
        while (g.R.size()) pushR(g.popR());
        while (topL() > topR()) {
            SLOPE l = popL(), r = popR();
            f0 -= (l - r);
            pushL(r), pushR(l);
        }
        return *this;
    }
};
