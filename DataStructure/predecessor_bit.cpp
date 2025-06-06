//
// BIT 木上の二分探索を用いた multiset の実装
//   ・クエリ先読みができる場合に有効
//   ・正の整数 M (10^7 程度) に対して、M 以下の非負整数を管理する multiset
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
//   ARC 197 C - Removal of Multiples
//     https://atcoder.jp/contests/abc356/tasks/abc356_f
//
//   ABC 229 G - Longest Y
//     https://atcoder.jp/contests/abc229/tasks/abc229_g
//


#include <bits/stdc++.h>
using namespace std;


// multiset by BIT
// manage integers x (0 <= x < lim), where lim <= 10^7
// Abel: type of the number of inserted values (not the tyep of the value)
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
    void insert(int x, Abel num = 1) { add(x, num); }
    void erase(int x, Abel num = 1) { add(x, -min(num, count(x))); }
    void clear() { dat.assign(lim, IDENTITY); }
    Abel count(int x) const { return sum(x, x + 1); }
    Abel count(int l, int r) const { return sum(l, r); }
    Abel size() const { return sum(lim); }
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
    int operator [] (int k) const { return get(k); }
    
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
        return prev(lim);
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

// ARC 197 C - Removal of Multiples
void ARC_197_C() {
    const int MAX = 3000000;
    FastMultiSetByBIT<int> bit(MAX);
    for (int v = 1; v < MAX; v++) bit.insert(v);

    int Q, A, B;;
    cin >> Q;
    for (int i = 0; i < Q; i++) {
        cin >> A >> B, B--;
        if (bit.count(A)) {
            for (int v = A; v < MAX; v += A) bit.erase(v);
        }
        cout << bit[B] << '\n';
    }
}

// ABC 229 G - Longest Y
void ABC_229_G() {
    using ll = long long;
    string S;
    ll K;
    cin >> S >> K;
    vector<ll> vs;
    for (int i = 0; i < S.size(); i++) if (S[i] == 'Y') vs.push_back(i - (int)vs.size());
    if (vs.size() <= 1) {
        cout << vs.size() << endl;
        return;
    }

    const int MAX = 210000;
    FastMultiSetByBIT<ll> se(MAX);
    auto solve = [&](ll x) -> ll {
        ll pre_sum = 0, suf_sum = 0, pre_num = x/2, suf_num = x - x/2;
        se.clear();
        for (int i = 0; i < x/2; i++) {
            se.insert(vs[i]);
            pre_sum += vs[i];
        }
        for (int i = x/2; i < x; i++) {
            se.insert(vs[i]);
            suf_sum += vs[i];
        }
        ll med = se[x/2];
        ll sum = (med * pre_num - pre_sum) + (suf_sum - med * suf_num);
        ll res = sum;
        for (int i = x; i < vs.size(); i++) {
            ll pre_before = se[0], pre_after = se[x/2];
            ll suf_before = se[x/2];
            se.insert(vs[i]), se.erase(vs[i-x]);
            ll suf_after = se[x-1];
            pre_sum -= pre_before, pre_sum += pre_after;
            suf_sum -= suf_before, suf_sum += suf_after;
            med = se[x/2];
            sum = (med * pre_num - pre_sum) + (suf_sum - med * suf_num);
            res = min(res, sum);
        }
        return res;
    };
    ll low = -1, high = (ll)vs.size() + 1;
    while (high - low > 1) {
        ll x = (high + low) / 2;
        if (solve(x) <= K) low = x;
        else high = x;
    }
    cout << low << endl;
}


int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    
    //Yosupo_Predecessor_Problem();
    //ARC_033_C();
    //ARC_197_C();
    ABC_229_G();
}