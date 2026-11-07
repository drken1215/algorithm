#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// find min_j f(i, j) for all i, by Monotone Minima, O(H + W log H)
// f(i, j) must be monotone (argmin is not decreasing)
// slide access OK
template<class VAL, class FUNC> void MonotoneMinimaRec
(int HL, int HR, int WL, int WR, const FUNC &f, vector<pair<VAL, int>> &res) {
    if (HR - HL <= 0) return;
    int HM = (HL + HR) / 2;
    res[HM].second = WL;
    for (int i = WL; i < WR; i++) {
        VAL val = f(HM, i);
        if (res[HM].first > val) res[HM] = make_pair(val, i);
    }
    MonotoneMinimaRec(HL, HM, WL, res[HM].second + 1, f, res);
    MonotoneMinimaRec(HM + 1, HR, res[HM].second, WR, f, res);
}
template<class VAL, class FUNC> vector<pair<VAL, int>> MonotoneMinima(int H, int W, const FUNC &f) {
    vector<pair<VAL, int>> res(H, make_pair(numeric_limits<VAL>::max() / 2, -1));
    MonotoneMinimaRec(0, H, 0, W, f, res);
    return res;
}
