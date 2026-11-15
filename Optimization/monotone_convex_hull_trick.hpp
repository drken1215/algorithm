#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// Convex Hull Trick
/*
    追加クエリの直線の傾きが単調減少の場合
 
    - insert(a, b): add y = ax + b, O(1)
    - query(x): min_i{a[i]*x + b}, O(log N)
    - query_monotone(x): min_i{a[i]*x + b}, O(1), x も単調増加であることを仮定
*/
template<class T> struct MonotoneConvexHullTrick {
    using Line = pair<T, T>;
    deque<Line> lines;
    T INF;

    MonotoneConvexHullTrick(T inf = numeric_limits<T>::max() / 2) : INF(inf) {}
    bool check(const Line &a, const Line &b, const Line &c) {
        return (b.first - a.first) * (c.second - b.second)
        >= (b.second - a.second) * (c.first - b.first);
    }
    
    T get(int k, T x) {
        if (lines.empty() || k >= lines.size()) return INF;
        assert(k >= 0 && k < lines.size());
        return lines[k].first * x + lines[k].second;
    }
    
    void insert(T a, T b) {
        Line l(a, b);
        while (lines.size() >= 2 && check(lines[(int)lines.size()-2], lines[(int)lines.size()-1], l)) {
            lines.pop_back();  // 不要な直線を捨てる
        }
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
        if (lines.empty()) return INF;
        return lines[0].first * x + lines[0].second;
    }
};
