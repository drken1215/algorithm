#define PROBLEM "https://onlinejudge.u-aizu.ac.jp/problems/DPL_3_C"

#include <bits/stdc++.h>
using namespace std;

#include "histogram.hpp"

int main() {
    int N;
    cin >> N;
    vector<long long> H(N);
    for (int i = 0; i < N; ++i) cin >> H[i];
    cout << max_area_in_histogram(H) << endl;
}
