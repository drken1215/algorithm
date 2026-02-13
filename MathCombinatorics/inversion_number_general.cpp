//
// 多重集合として一致する 2　系列 A, B 間の転倒数を求める
//   ・同じ数値は、A, B それぞれ左から順に対応をとればよい
//
// verified
//   AOJ Course ALDS1_5_D Recursion / Divide and Conquer - The Number of Inversions
//     http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ALDS1_5_D&lang=jp
//


#include <bits/stdc++.h>
using namespace std;


// Fast MultiSet By BIT
template<class Abel = int> struct FastMultiSetByBIT {
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

// mapping[i]: A[i] に対応する B の要素の index
template<class T> vector<long long> find_mapping(vector<T> A, vector<T> B) {
    // 多重集合として等しいことを保証して、座標圧縮する
    int N = (int)A.size();
    auto A2 = A, B2 = B;
    sort(A2.begin(), A2.end()), sort(B2.begin(), B2.end());
    assert(A2 == B2);
    A2.erase(unique(A2.begin(), A2.end()), A2.end());
    for (int i = 0; i < N; i++) {
        A[i] = lower_bound(A2.begin(), A2.end(), A[i]) - A2.begin();
        B[i] = lower_bound(A2.begin(), A2.end(), B[i]) - A2.begin();
    }

    // B の各値ごとに index を求める
    vector<vector<long long>> pb(N);
    for (int i = 0; i < N; i++) pb[B[i]].emplace_back(i);

    // A[i] が B で何番目なのかを求める
    vector<long long> res(N), iter(N, 0);
    for (int i = 0; i < N; i++) res[i] = pb[A[i]][iter[A[i]]++];
    return res;
}

// A の転倒数
template<class T> T inversion_number(vector<T> A) {
    int N = (int)A.size();
    auto A2 = A;
    sort(A2.begin(), A2.end());
    A2.erase(unique(A2.begin(), A2.end()), A2.end());
    for (int i = 0; i < N; i++) A[i] = lower_bound(A2.begin(), A2.end(), A[i]) - A2.begin();

    T res = 0;
    FastMultiSetByBIT<T> S(N);
    for (int i = 0; i < N; i++) {
        res += S.count(A[i] + 1, N);
        S.add(A[i], 1);
    }
    return res;
}

// A, B の間の転倒数
template<class T> T inversion_number(vector<T> A, vector<T> B) {
    auto mapping = find_mapping(A, B);
    return inversion_number(mapping);
}


//------------------------------//
// Examples
//------------------------------//

int main() {
    int n; cin >> n;
    vector<int> a(n); for (int i = 0; i < n; ++i) cin >> a[i];
    cout << inversion_number(a) << endl;
}
