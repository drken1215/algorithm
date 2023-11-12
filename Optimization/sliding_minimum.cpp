//
// スライド最小値
//
// cf.
//
//
// verified:
//   AOJ 0613 財宝 (JOI 2014 予選 F)
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=0613
//

#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
using namespace std;
using pll = pair<long long, long long>;


// slide min
template<class VAL, class TIME> struct SlideMin {
    const VAL INF = 1LL<<60; // to be set appropriately
    const TIME nul = -1; // to be set appropriately
        
    deque<pair<VAL, TIME> > deq;
    SlideMin()  { }

    // get minimum
    pair<VAL, TIME> get() {
        if (deq.empty()) return make_pair(INF, nul);
        return deq.front();
    }

    // push (v, t), "t is non-decreasing" is necessary
    void push(VAL v, TIME t) {
        while (!deq.empty() && deq.back().first >= v) deq.pop_back();
        deq.emplace_back(v, t);
    }

    // pop "< t", "t it non-decreasing" is necessary
    void pop(TIME t) {
        while (!deq.empty() && deq.front().second < t) deq.pop_front();
    }
};



vector<long long> X, Y;
void rec(int i, int last, long long x, long long y, vector<pll> &v) {
    if (i == last) {
        v.push_back(pll(x, y));
        return;
    }
    rec(i + 1, last, x, y, v);
    rec(i + 1, last, x + X[i], y - Y[i], v);
    rec(i + 1, last, x - X[i], y + Y[i], v);
}

int main() {
    int N; long long D; cin >> N >> D;
    X.resize(N); Y.resize(N);
    for (int i = 0; i < N; ++i) cin >> X[i] >> Y[i];
    vector<pll> L, R;
    rec(0, N/2, 0, 0, L);
    rec(N/2, N, 0, 0, R);
    sort(L.begin(), L.end(), greater<pll>());
    sort(R.begin(), R.end());

    long long res = 0;
    SlideMin<long long, long long> sm;
    int pos = 0;
    for (auto it : L) {
        // push
        while (pos < R.size() && R[pos].first <= -it.first + D) {
            sm.push(-R[pos].second, R[pos].first);
            ++pos;
        }
        
        // pop
        sm.pop(-it.first - D);

        // get
        long long minv = -sm.get().first;
        res = max(res, it.second + minv);
    }
    cout << res << endl;
}
