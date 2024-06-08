//
// BIT 木上の二分探索を用いた Multiset の実装
//
// verified:
//   Yosupo Library Checker - Predecessor Problem
//     https://judge.yosupo.jp/problem/predecessor_problem
//
//   ARC 033 C - データ構造
//     https://beta.atcoder.jp/contests/arc033/tasks/arc033_3
//
//   ABC 356 F - Distance Component Size Query
//     https://atcoder.jp/contests/abc356/tasks/abc356_f
//


#include <bits/stdc++.h>
using namespace std;


template<class Abel> struct FastMultiSetByBIT {
    int topbit(int x) const { return (x == 0 ? -1 : 31 - __builtin_clz(x)); }
    int lowbit(int x) const { return (x == 0 ? -1 : __builtin_ctz(x)); }
    int lim;
    Abel IDENTITY;
    vector<Abel> dat;
    
    // [0, n)
    FastMultiSetByBIT(int n, Abel identity = 0)
    : lim(n), IDENTITY(identity), dat(n, identity) { }
    void init(int n, Abel identity = 0) {
        lim = n;
        IDENTITY = identity;
        dat.assign(n, IDENTITY);
    }
    
    // p is 0-indexed
    void add(int p, Abel x) {
        if (p < 0) p = 0;
        for (int i = p; i < (int)dat.size(); i |= i + 1)
            dat[i] = dat[i] + x;
    }
    
    // [0, p), p is 0-indexed
    Abel sum(int p) const {
        if (p > lim) p = lim;
        Abel res = IDENTITY;
        for (int i = p - 1; i >= 0; i = (i & (i + 1)) - 1)
            res = res + dat[i];
        return res;
    }
    
    // [l, r), l and r are 0-indexed
    Abel sum(int l, int r) const {
        return sum(r) - sum(l);
    }
    
    // insert, erase, count, min, max
    void insert(int x) { add(x, 1); }
    void erase(int x) { if (count(x)) add(x, -1); }
    Abel count(int x) const { return sum(x, x + 1); }
    Abel count(int l, int r) const { return sum(l, r); }
    Abel size() const { return sum(lim); }
    bool operator [] (int x) const { return count(x); }
    int get_min() const { return next(); }
    int get_max() const { return prev(); }

    // get max r s.t. check(sum(l, r)) = True (0-indexed), O(log N)
    // check(IDENTITY) must be True
    int max_right(const function<bool(Abel)> check, int l = 0) const {
        if (l >= lim) return lim;
        assert(check(IDENTITY));
        Abel s = IDENTITY;
        int k = 0;
        while (true) {
            if (l % 2 == 1) s = s - dat[l - 1], --l;
            if (l <= 0) {
                k = topbit(lim) + 1;
                break;
            }
            k = lowbit(l) - 1;
            if (l + (1 << k) > lim) break;
            if (!check(s + dat[l + (1 << k) - 1])) break;
            s = s - dat[l - 1];
            l -= l & -l;
        }
        while (k) {
            --k;
            if (l + (1 << k) - 1 < lim) {
                Abel ns = s + dat[l + (1 << k) - 1];
                if (check(ns)) {
                    l += (1 << k);
                    s = ns;
                }
            }
        }
        return l;
    }
    
    // get min l s.t. check(sum(l, r))  = True (0-indexed), O(log N)
    // check(IDENTITY) must be True
    int min_left(const function<bool(Abel)> check, int r = -1) const {
        if (r == -1) r = lim;
        if (r <= 0) return 0;
        assert(check(IDENTITY));
        Abel s = IDENTITY;
        int k = 0;
        while (r > 0 && check(s)) {
            s = s + dat[r - 1];
            k = lowbit(r);
            r -= r & -r;
        }
        if (check(s)) return 0;
        while (k) {
            --k;
            Abel ns = s - dat[r + (1 << k) - 1];
            if (!check(ns)) {
                r += (1 << k);
                s = ns;
            }
        }
        return r + 1;
    }
              
    // k-th number that is not less than l (k is 0-indexed)
    int get(Abel k, int l = 0) const {
        return max_right([&](Abel x) { return x <= k; }, l);
    }
    
    // next (including x)
    int next(int l = 0) const {
        if (l < 0) l = 0;
        if (l > lim) l = lim;
        return max_right([&](Abel x) { return x <= 0; }, l);
    }
    
    // prev (including x)
    int prev(int r) const {
        if (r > lim) r = lim;
        return min_left([&](Abel x) { return x <= 0; }, r + 1) - 1;
    }
    int prev() const {
        prev(lim);
    }
    
    // debug
    friend ostream& operator << (ostream &s, const FastMultiSetByBIT &fs) {
        for (int x = fs.get_min(); x < fs.lim; x = fs.next(x + 1)) {
            s << x << " ";
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
    FastMultiSetByBIT<int> fs(N);
    for (int i = 0; i < N; ++i) if (T[i] == '1') fs.insert(i);
    while (Q--) {
        int t, k;
        cin >> t >> k;
        if (t == 0) {
            if (!fs.count(k)) fs.insert(k);
        } else if (t == 1) {
            if (fs.count(k)) fs.erase(k);
        } else if (t == 2) cout << (fs.count(k) ? 1 : 0) << '\n';
        else if (t == 3) {
            int res = fs.next(k);
            cout << (res < N ? res : -1) << '\n';
        }
        else cout << fs.prev(k) << '\n';
    }
}

// ARC 033 C - データ構造
void ARC_033_C() {
    FastMultiSetByBIT<int> fs(200000);
    int Q;
    cin >> Q;
    for (int i = 0; i < Q; ++i) {
        int T, X;
        cin >> T >> X;
        if (T == 1) fs.insert(X); // insert X
        else {
            int val = fs.get(X - 1); // X is 1-indexed
            cout << val << endl;
            fs.erase(val); // delete val
        }
    }
}


int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    
    //Yosupo_Predecessor_Problem();
    ARC_033_C();
}
