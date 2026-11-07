#pragma once
#include <bits/stdc++.h>
#include <cassert>
using namespace std;

// disjoint sparse table
template<class SemiGroup> struct DisjointSparseTable {
    using Func = function< SemiGroup(SemiGroup, SemiGroup) >;
    const Func F;
    vector<vector<SemiGroup> > dat;
    vector<int> height;
    
    DisjointSparseTable(const Func &f) : F(f) { }
    DisjointSparseTable(const Func &f, const vector<SemiGroup> &vec) : F(f) { init(vec); }
    void init(const vector<SemiGroup> &vec) {
        int n = (int)vec.size(), h = 1;
        while ((1<<h) <= n) ++h;
        dat.assign(h, vector<SemiGroup>(n));
        height.assign((1<<h), 0);
        for (int i = 2; i < (1<<h); i++) height[i] = height[i>>1]+1;
        for (int i = 0; i < n; ++i) dat[0][i] = vec[i];
        for (int i = 1; i < h; ++i) {
            int s = (1<<i);
            for (int j = 0; j < n; j += (s << 1)) {
                int t = min(j+s, n);
                dat[i][t-1] = vec[t-1];
                for (int k = t-2; k >= j; --k) dat[i][k] = F(vec[k], dat[i][k+1]);
                if (n <= t) break;
                dat[i][t] = vec[t];
                for (int k = t+1; k < min(t+s, n); ++k) dat[i][k] = F(dat[i][k-1], vec[k]);
            }
        }
    }
    SemiGroup get(int a, int b) {
        if (a >= --b) return dat[0][a];
        return F(dat[height[a^b]][a], dat[height[a^b]][b]);
    }
};
