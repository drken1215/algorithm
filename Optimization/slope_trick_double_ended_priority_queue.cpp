//
// Slope Trick (両端 priority queue で実装, 通常の priority queue よりもメモリを食うので注意)
//
// verified
//   第2回 ドワンゴからの挑戦状 予選 E - 花火 (for 累積min, +|x-a|)
//     https://atcoder.jp/contests/dwango2016-prelims/tasks/dwango2016qual_e
//
//   AtCoder ABC 217 H - Snuketoon (for slide(-dt, dt), +max(0, a-x) など)
//     https://atcoder.jp/contests/abc217/tasks/abc217_h
//
//   KUPC 2016 H - 壁壁壁壁壁壁壁 (for slide(-INF, A[i]-B[i]), eval(x) など)
//     https://atcoder.jp/contests/kupc2016/tasks/kupc2016_h
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


// Double Ended Priority Queue
template<class T> struct DoubleEndedPriorityQueue {
    // removable heap
    template<class QUETYPE> struct RemovablePriorityQueue {
        using VALTYPE = typename QUETYPE::value_type;
        QUETYPE que, delay;
        
        // constructor
        RemovablePriorityQueue() {}
        RemovablePriorityQueue(const RemovablePriorityQueue&) = default;
        RemovablePriorityQueue& operator = (const RemovablePriorityQueue&) = default;
        
        // getter
        constexpr int size() const { return (int)que.size() - (int)delay.size(); }
        constexpr bool empty() const { return size() == 0; }

        // push(x), remove(x)
        constexpr void push(VALTYPE x) { que.push(x); }
        constexpr void remove(VALTYPE x) { delay.push(x); }
        
        // pop min/max value
        constexpr VALTYPE pop() {
            T res = get();
            que.pop();
            return res;
        }
        
        // get min/max value (not pop)
        constexpr VALTYPE get() {
            assert(!que.empty());
            while (!delay.empty() && que.top() == delay.top()) {
                que.pop();
                delay.pop();
            }
            assert(!que.empty());
            return que.top();
        }
    };
    
    // inner data
    RemovablePriorityQueue<priority_queue<T, vector<T>, greater<T>>> min_que;
    RemovablePriorityQueue<priority_queue<T>> max_que;
    
    // constructor
    DoubleEndedPriorityQueue() {}
    DoubleEndedPriorityQueue(const DoubleEndedPriorityQueue&) = default;
    DoubleEndedPriorityQueue& operator = (const DoubleEndedPriorityQueue&) = default;
    
    // getter
    constexpr int size() const {
        return (int)min_que.size();
    }
    constexpr bool empty() const {
        return size() == 0;
    }
    
    // push(x), remove(x)
    constexpr void push(T x) {
        min_que.push(x);
        max_que.push(x);
    }
    constexpr void remove(T x) {
        min_que.remove(x);
        max_que.remove(x);
    }
    
    // get min, pop min
    constexpr T get_min() {
        return min_que.get();
    }
    constexpr T pop_min() {
        T x = min_que.pop();
        max_que.remove(x);
        return x;
    }
    
    // get max, pop max
    constexpr T get_max() {
        return max_que.get();
    }
    constexpr T pop_max() {
        T x = max_que.pop();
        min_que.remove(x);
        return x;
    }
};

// Slope Trick
/*
    f(x): 区分線形凸関数
    ・min_f: f(x) の最小値
    ・offset_l, offset_r: L, R 全体に対して加算する値
    ・L, R: f(x) の傾きが 1 変化する x 座標の多重集合 (min_f の左側と右側)
*/
template<class COORD> struct SlopeTrickByEDPQ {
    using POINT = pair<COORD, COORD>;
    using LINE = vector<pair<COORD, COORD>>;

    // inner data
    COORD min_f, offsetL, offsetR, INF;
    DoubleEndedPriorityQueue<COORD> L, R;

    // constructors
    SlopeTrickByEDPQ(COORD inf = numeric_limits<COORD>::max() / 2)
        : min_f(0), offsetL(0), offsetR(0), INF(inf) {
        assert(inf > 0);
    }
    SlopeTrickByEDPQ(const SlopeTrickByEDPQ&) = default;
    SlopeTrickByEDPQ& operator = (const SlopeTrickByEDPQ&) = default;

    // basic operations
    constexpr int sizeL() const { return (int)L.size(); }
    constexpr int sizeR() const { return (int)R.size(); }
    constexpr int size() const { return sizeL() + sizeR(); }
    constexpr void pushL(const COORD &v) { L.push(v - offsetL); }
    constexpr void pushR(const COORD &v) { R.push(v - offsetR); }
    constexpr COORD topL() { return L.empty() ? -INF : L.get_max() + offsetL; }
    constexpr COORD topR() { return R.empty() ? INF : R.get_min() + offsetR; }
    constexpr COORD popL() {
        auto res = topL();
        if (!L.empty()) L.pop_max();
        return res;
    }
    constexpr COORD popR() {
        auto res = topR();
        if (!R.empty()) R.pop_min();
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
            auto t = L2.get_max() + offsetL;
            L2.pop_max();
            res += max(COORD(0), t - x);
        }
        while (!R2.empty()) {
            auto t = R2.get_min() + offsetR;
            R2.pop_min();
            res += max(COORD(0), x - t);
        }
        return res + min_f;
    }
    constexpr friend ostream &operator << (ostream &os, SlopeTrickByEDPQ st) {
        auto [lineL, lineR] = st.get_lines();
        os << endl << "left: ";
        for (auto [x, y] : lineL) os << "(" << x << ", " << y << ") ";
        os << endl << "right: ";
        for (auto [x, y] : lineR) os << "(" << x << ", " << y << ") ";
        return os << endl;
    }
    
    // f(x) += b, O(1)
    SlopeTrickByEDPQ &add_const(const COORD &b) {
        min_f += b;
        return *this;
    }

    // f(x) += max(0, x - a), O(log N)
    SlopeTrickByEDPQ &add_relu(const COORD &a) {
        if (topL() <= a) pushR(a);
        else {
            min_f += max(COORD(0), topL() - a);
            pushL(a), pushR(popL());
        }
        return *this;
    }

    // f(x) += max(0, a - x), O(log N)
    SlopeTrickByEDPQ &add_irelu(const COORD &a) {
        if (topR() >= a) pushL(a);
        else {
            min_f += max(COORD(0), a - topR());
            pushR(a), pushL(popR());
        }
        return *this;
    }

    // f(x) += |x - a|, O(log N)
    SlopeTrickByEDPQ &add_abs(const COORD &a) {
        add_relu(a), add_irelu(a);
        return *this;
    }

    // f(x) <- g(x) = min_{y <= x} f(y), (\_/ -> \__), O(1)
    SlopeTrickByEDPQ &clear_right() {
        R = DoubleEndedPriorityQueue<COORD>();
        return *this;
    }

    // f(x) <- g(x) = min_{y >= x} f(y), (\_/ -> __/), O(1)
    SlopeTrickByEDPQ &clear_left() {
        L = DoubleEndedPriorityQueue<COORD>();
        return *this;
    }

    // f(x) <- 0, O(1)
    SlopeTrickByEDPQ &clear() {
        *this = SlopeTrickByEDPQ();
        return *this;
    }

    // f(x) <- g(x) = f(x - a), O(1)
    SlopeTrickByEDPQ &slide(const COORD &a) {
        offsetL += a, offsetR += a;
        return *this;
    }

    // f(x) <- g(x) = min_{a <= y <= b} f(x - y) = min_{x-b <= y <= x-a} f(y), O(1)
    SlopeTrickByEDPQ &slide(const COORD &a, const COORD &b) {
        assert(a <= b);
        if (a <= -INF) clear_left();
        else offsetL += a;
        if (b >= INF) clear_right();
        else offsetR += b;
        return *this;
    }

    // f(x) <- g(x) = min_{0 <= y <= a} f(x - y) = min_{x-a <= y <= x} f(y), O(1)
    SlopeTrickByEDPQ &slide_right_curve_to_right(const COORD &a) {
        assert(a >= 0);
        slide(0, a);
        return *this;
    }

    // f(x) <- g(x) = min_{-a <= y <= 0} f(x - y) = min_{x <= y <= x+a} f(y),  O(1)
    SlopeTrickByEDPQ &slide_left_curve_to_left(const COORD &a) {
        assert(a >= 0);
        slide(-a, 0);
        return *this;
    }

    // f(x) += g(x), O((log N)^2)
    // attention: this function is destructive
    SlopeTrickByEDPQ &add(SlopeTrickByEDPQ &g) {
        if (size() < g.size()) {
            swap(min_f, g.min_f);
            swap(offsetL, g.offsetL), swap(offsetR, g.offsetR);
            swap(L, g.L), swap(R, g.R);
        }
        min_f += g.min_f;
        while (g.L.size()) add_irelu(g.popL());
        while (g.R.size()) add_relu(g.popR());
        g.clear();
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
    SlopeTrickByEDPQ<long long> dp;
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
    SlopeTrickByEDPQ<long long> dp;
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
    SlopeTrickByEDPQ<long long> dp;
    for (int i = 0; i < N * 2; i++) dp.add_abs(0);
    for (int i = 0; i < N; i++) {
        dp.slide(-dp.INF, A[i]-B[i]);
        dp.add_abs(0);
    }
    long long res = dp.eval(0);
    cout << res << endl;
}


int main() {
    //DWANGO_2nd_prelims_E();
    //ABC_217_H();
    KUPC_2016_H();
}