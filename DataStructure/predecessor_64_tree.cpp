//
// 64-ary tree
//
// verified:
//   Yosupo Library Checker - Predecessor Problem
//     https://judge.yosupo.jp/problem/predecessor_problem
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#include <bits/stdc++.h>
using namespace std;


// 64-ary tree
// manage integers x (0 <= x < lim), where lim <= 10^7
struct FastSet {
    using u32 = unsigned int;
    using u64 = unsigned long long ;
    int lowbit(u64 x) const { return (x == 0 ? -1 : __builtin_ctzll(x)); }
    int topbit(u64 x) const { return (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
    static constexpr u32 BASE = 64;
    int lim, log;
    vector<vector<u64>> seg;
    
    // constructor
    FastSet() {}
    FastSet(int n) { build(n); }
    template<class F> FastSet(int n, F isin) { build(n, isin); }
    void build(int n) {
        lim = n;
        seg.clear();
        do {
            seg.push_back(vector<u64>((n + BASE - 1) / BASE));
            n = (n + BASE - 1) / BASE;
        } while (n > 1);
        log = (int)seg.size();
    }
    template<class F> void build(int n, F isin) {
        build(n);
        for (int x = 0; x < n; ++x) {
            seg[0][x / BASE] |= u64(isin(x)) << (x % BASE);
        }
        for (int h = 0; h < log - 1; ++h) {
            for (int x = 0; x < (int)seg[h].size(); ++x) {
                seg[h + 1][x / BASE] |= u64(bool(seg[h][x])) << (x % BASE);
            }
        }
    }
    
    // insert / erase / exist
    void insert(int x) {
        for (int h = 0; h < log; ++h) {
            seg[h][x / BASE] |= u64(1) << (x % BASE), x /= BASE;
        }
    }
    void erase(int x) {
        u64 y = 0;
        for (int h = 0; h < log; ++h) {
            seg[h][x / BASE] &= ~(u64(1) << (x % BASE));
            seg[h][x / BASE] |= y << (x % BASE);
            y = bool(seg[h][x / BASE]);
            x /= BASE;
        }
    }
    bool isin(int x) const {
        return seg[0][x / BASE] >> (x % BASE) & 1;
    }
    bool operator [] (int x) const { return isin(x); }
    
    
    // next / prev
    int next(int x) const {
        assert(x <= lim);
        if (x < 0) x = 0;
        for (int h = 0; h < log; ++h) {
            if (x / BASE == seg[h].size()) break;
            u64 div = seg[h][x / BASE] >> (x % BASE);
            if (!div) {
                x = x / BASE + 1;
                continue;
            }
            x += lowbit(div);
            for (int g = h - 1; g >= 0; --g) {
                x *= BASE;
                x += lowbit(seg[g][x / BASE]);
            }
            return x;
        }
        return -1;  // not exist
    }
    int prev(int x) const {
        assert(x >= -1);
        if (x >= lim) x = lim - 1;
        for (int h = 0; h < log; ++h) {
            if (x == -1) break;
            u64 div = seg[h][x / BASE] << (63 - x % BASE);
            if (!div) {
                x = x / BASE - 1;
                continue;
            }
            x -= __builtin_clzll(div);
            for (int g = h - 1; g >= 0; --g) {
                x *= BASE;
                x += topbit(seg[g][x / BASE]);
            }
            return x;
        }
        return -1;  // not exist
    }
    
    // debug
    friend ostream& operator << (ostream &s, const FastSet &fs) {
        int x = fs.next(0);
        while (x != -1) {
            s << x << " ";
            x = fs.next(x + 1);
        }
        return s;
    }
};



//------------------------------//
// Examples
//------------------------------//

// Yosupo Library Checker - Predecessor Problem
void Yosupo_Predecessor_Problem() {
    int N, Q;
    string T;
    cin >> N >> Q >> T;
    FastSet fs(N, [&](int x){ return T[x] == '1'; });
    while (Q--) {
        int t, k;
        cin >> t >> k;
        if (t == 0) fs.insert(k);
        else if (t == 1) fs.erase(k);
        else if (t == 2) cout << fs.isin(k) << '\n';
        else if (t == 3) cout << fs.next(k) << '\n';
        else cout << fs.prev(k) << '\n';
        
        // cout << fs << '\n';
    }
}


int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    
    Yosupo_Predecessor_Problem();
}

