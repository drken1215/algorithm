//
// Slope Trick
//
// verified
//   第2回 ドワンゴからの挑戦状 予選 E - 花火
//     https://atcoder.jp/contests/dwango2016-prelims/tasks/dwango2016qual_e
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

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


// Double Ended Priority Queue
template<class T> struct DoubleEndedPriorityQueue {
    // removable heap
    template<class QUETYPE> struct RemovablePriorityQueue {
        using VALTYPE = typename QUETYPE::value_type;
        QUETYPE que, delay;
        
        // constructor
        RemovablePriorityQueue() {}
        
        // getter
        int size() { return (int)que.size() - (int)delay.size(); }
        bool empty() { return size() == 0; }

        // push(x), remove(x)
        void push(VALTYPE x) { que.push(x); }
        void remove(VALTYPE x) { delay.push(x); }
        
        // pop min/max value
        VALTYPE pop() {
            T res = get();
            que.pop();
            return res;
        }
        
        // get min/max value (not pop)
        VALTYPE get() {
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
    
    // getter
    int size() {
        return (int)min_que.size();
    }
    bool empty() {
        return size() == 0;
    }
    
    // push(x), remove(x)
    void push(T x) {
        min_que.push(x);
        max_que.push(x);
    }
    void remove(T x) {
        min_que.remove(x);
        max_que.remove(x);
    }
    
    // get min, pop min
    T get_min() {
        return min_que.get();
    }
    T pop_min() {
        T x = min_que.pop();
        max_que.remove(x);
        return x;
    }
    
    // get max, pop max
    T get_max() {
        return max_que.get();
    }
    T pop_max() {
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
template<class COORD> struct SlopeTrick {
    // inner data
    COORD min_f, offsetL, offsetR, INF;
    DoubleEndedPriorityQueue<COORD> L, R;

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
    constexpr void pushL(const COORD &v) const { L.push(v - offsetL); }
    constexpr void pushR(const COORD &v) const { R.push(v - offsetR); }
    constexpr COORD topL() const { return L.empty() ? -INF : L.get_max() - offsetL; }
    constexpr COORD topR() const { return R.empty() ? -INF : R.get_min() - offsetR; }
    constexpr COORD popL() const {
        auto res = topL();
        if (!L.empty()) L.pop_max();
        return res;
    }
    constexpr COORD popR() const {
        auto res = topR();
        if (!R.empty()) R.pop_min();
        return res;
    }
    
    // f(x) += b, O(1)
    SlopeTrick &add_const(const COORD &b) {
        min_f += b;
        return *this;
    }

    // f(x) += max(0, x - a), O(log N)
    SlopeTrick &add_relu(const COORD &a) {
        min_f += max(COORD(0), topL() - a);
        pushL(a), pushR(popL());
        return *this;
    }

    // f(x) += max(0, a - x), O(log N)
    SlopeTrick &add_irelu(const COORD &a) {
        min_f += max(COORD(0), a - topR());
        pushR(a), pushL(popR());
        return *this;
    }

    // f(x) += |x - a|, O(log N)
    SlopeTrick &add_abs(const COORD &a) {
        add_relu(a), add_irelu(a);
        return *this;
    }

    // f(x) <- g(x) = min_{y <= x} f(y), (\_/ -> \__), O(1)
    SlopeTrick &clear_right() {
        R.clear();
        return *this;
    }

    // f(x) <- g(x) = min_{y >= x} f(y), (\_/ -> __/), O(1)
    SlopeTrick &clear_left() {
        L.clear();
        return *this;
    }

    // f(x) <- g(x) = f(x - a), O(1)
    SlopeTrick &slide(const COORD &a) {
        offsetL += a, offsetR += a;
        return *this;
    }

    // f(x) <- g(x) = min_{a <= y <= b} f(x - y)
    SlopeTrick &slide(const COORD &a, const COORD &b) {
        assert(a <= b);
        offsetL += a, offsetR += b;
        return *this;
    }

    // f(x) <- g(x) = min_{0 <= y <= a} f(x - y)
    SlopeTrick &slide_right_curve_to_right(const COORD &a) {
        assert(a >= 0);
        offsetR += a;
        return *this;
    }

    // f(x) <- g(x) = min_{0 <= y <= a} f(x + y)
    SlopeTrick &slide_left_curve_to_left(const COORD &a) {
        assert(a >= 0);
        offsetL += a;
        return *this;
    }


};


//------------------------------//
// Solver
//------------------------------//

// 第2回 ドワンゴからの挑戦状 予選 E - 花火
void DWANGO_2nd_prelims_E() {

}


int main() {
    DWANGO_2nd_prelims_E();
}