#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// min value of H within [left[i], right[i]) <= H[i]
pair<vector<int>, vector<int>> solve_left_right(const vector<long long> &H) {
    int N = (int)H.size();
    vector<int> left(N, 0), right(N, N);
    
    // left
    stack<pair<long long, int>> stack_left;
    for(int i = 0; i < N; ++i) {
        while (!stack_left.empty() && H[i] <= stack_left.top().first)
            stack_left.pop();
        if (!stack_left.empty()) left[i] = stack_left.top().second + 1;
        stack_left.push({H[i], i});
    }
    
    // right
    stack<pair<long long, int>> stack_right;
    for(int i = N-1; i >= 0; --i) {
        while (!stack_right.empty() && H[i] <= stack_right.top().first)
            stack_right.pop();
        if (!stack_right.empty()) right[i] = stack_right.top().second;
        stack_right.push({H[i], i});
    }
    return {left, right};
}

// max area of rectangle in histogram
long long max_area_in_histogram(const vector<long long> &H) {
    auto [left, right] = solve_left_right(H);
    long long res = 0;
    for (int i = 0; i < H.size(); ++i) {
        res = max(res, H[i] * (right[i] - left[i]));
    }
    return res;
}
