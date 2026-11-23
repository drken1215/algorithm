#pragma once
#include <bits/stdc++.h>
#include <cassert>
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

    // basic operations
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

    // getter and debugger
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

    // f(x) <- g(x) = min_{y >= x} f(y), (\_/ -> __/), O(1)
    SlopeTrick &clear_left() {
        L = priority_queue<COORD>{};
        return *this;
    }

    // f(x) <- g(x) = min_{y <= x} f(y), (\_/ -> \__), O(1)
    SlopeTrick &clear_right() {
        R = priority_queue<COORD, vector<COORD>, greater<COORD>>{};
        return *this;
    }

    // f(x) <- 0, O(1)
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

    // f(x) += g(x), amotized O((log N)^2)
    // attention: this function is destructive to g
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

    // f(x) <- h(x) = min_{x = y + z, z >= 0} (f(y) + |z - a|), O(log N)
    SlopeTrick &min_plus_convolution_abs(const COORD &a) {
        assert(a >= 0);
        if (sizeL() == 0 || sizeR() == 0) return *this;
        auto l = topL(), r = topR();
        popL();
        pushL(l + a);
        R = priority_queue<COORD, vector<COORD>, greater<COORD>>{};
        pushR(r + a);
        return *this;
    }
};
