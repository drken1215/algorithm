//
// Convex Hull Trick
//
// verified
//   COLOCON 2018 Final C - スペースエクスプローラー高橋君
//     https://beta.atcoder.jp/contests/colopl2018-final-open/tasks/colopl2018_final_c
//


#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
using namespace std;


// Convex Hull Trick
/*
    追加クエリの直線の傾きが単調減少の場合
 
    - insert (a, b): add y = ax + b
    - query (x): min_i{a[i]*x + b}
*/
template<class T> struct CHT {
    using Line = pair<T, T>;
    deque<Line> lines;
    
    inline bool check(const Line &a, const Line &b, const Line &c) {
        return (b.first - a.first) * (c.second - b.second)
        >= (b.second - a.second) * (c.first - b.first);
    }
    
    inline T get(int k, T x) {
        return lines[k].first * x + lines[k].second;
    }
    
    void insert(T a, T b) {
        Line l(a, b);
        while (lines.size() >= 2
               && check(lines[(int)lines.size()-2], lines[(int)lines.size()-1], l))
            lines.pop_back();
        lines.push_back(l);
    }
    
    T query(T x) {
        int low = -1, high = (int)lines.size();
        while (high - low > 1) {
            int mid = (low + high) / 2;
            if (get(mid, x) >= get(mid + 1, x)) low = mid;
            else high = mid;
        }
        return get(high, x);
    }
    
    // クエリの単調性も成り立つ場合 (x が単調増加)
    T query_monotone(T x) {
        while (lines.size() >= 2 && get(0, x) >= get(1, x)) lines.pop_front();
        return lines[0].first * x + lines[0].second;
    }
};



int main() {
    long long N; cin >> N;
    vector<long long> a(N), res(N, 1LL<<60);
    for (int i = 0; i < N; ++i) cin >> a[i];
    CHT<long long> cht;
    for (long long i = 0; i < N; ++i) cht.insert(-2LL*i, a[i] + i*i);
    for (long long i = 0; i < N; ++i) res[i] = cht.query(i) + i*i;
    for (int i = 0; i < N; ++i) cout << res[i] << endl;
}
