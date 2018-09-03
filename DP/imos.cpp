//
// imos-algorithm
//
// c.f.
//   https://imoz.jp/algorithms/imos_method.html
//
// verified
//   AOJ Course DSL_5_A Cumulative Sum - The Maximum Number of Customers
//

#include <iostream>
#include <vector>
using namespace std;

typedef pair<int,int> Interval; // means the interval [first, second)

// T: max value of intervals
int imos(const vector<Interval> &intervals, int T) {
    vector<int> nums(T+1, 0);
    for (auto interval : intervals) {
        nums[interval.first]++;
        nums[interval.second]--;
    }
    for (int t = 0; t < T; ++t) {
        nums[t+1] += nums[t];
    }

    int res = 0;
    for (int t = 0; t <= T; ++t) res = max(res, nums[t]);

    return res;
}

int main() {
    int N, T; cin >> N >> T;
    vector<Interval> intervals(N);
    for (int i = 0; i < N; ++i) cin >> intervals[i].first >> intervals[i].second;
    cout << imos(intervals, T) << endl;
}  
