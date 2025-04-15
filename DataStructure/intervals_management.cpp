//
// 区間 (と値) を Set で管理する構造体
//
// verified
//   第六回 アルゴリズム実技検定 M - 等しい数
//     https://atcoder.jp/contests/past202104-open/tasks/past202104_m
//
//   RUPC 2018 G - Elevator
//     https://onlinejudge.u-aizu.ac.jp/problems/2880 
//


#include <bits/stdc++.h>
using namespace std;


// Interval Set
// T: type of range, VAL: data type
template<class T, class VAL = long long> struct IntervalSet {
    struct Node {
        T l, r;
        VAL val;
        Node(const T &l, const T &r, const VAL &val) : l(l), r(r), val(val) {}
        constexpr bool operator < (const Node &rhs) const {
            if (l != rhs.l) return l < rhs.l;
            else return r < rhs.r;
        }
    };

    // internal values
    const VAL identity;
    set<Node> S;

    // constructor
    IntervalSet(const VAL &identity = VAL()) : identity(identity) {}
    IntervalSet(const vector<VAL> &v, const VAL &identity = VAL()) : identity(identity) {
        vector<Node> vec;
        for (int l = 0; l < (int)v.size();) {
            int r = l;
            while (r < (int)v.size() && v[r] == v[l]) r++;
            vec.emplace_back(l, r, v[l]);
            l = r;
        }
        S = set<Node>(vec.begin(), vec.end());
    }

    // get the iterator of interval which contains p
    // not exist -> S.end()
    constexpr typename set<Node>::iterator get(const T &p) {
        auto it = S.upper_bound(Node(p, numeric_limits<T>::max(), 0));
        if (it == S.begin()) return S.end();
        it = prev(it);
        if (it->l <= p && p < it->r) return it;
        else return S.end();
    }

    // exist the interval which contains p: true
    constexpr bool contains(const T &p) {
        auto it = get(p);
        if (it != S.end()) return true;
        else return false;
    }

    // get the value of interval which contains p
    // not exist -> identity
    constexpr VAL get_val(const T &p) {
        auto it = get(p);
        if (it != S.end()) return it->val;
        else return identity;
    }
    VAL operator [] (const T &p) const {
        return get_val(p);
    }

    // get the leftist iterator of interval which contains value >= p
    constexpr typename set<Node>::iterator lower_bound(const T &p) {
        auto it = S.upper_bound(Node(p, numeric_limits<T>::max(), 0));
        if (it == S.begin()) return it;
        else return prev(it);
    }

    // is p, q in same interval?
    constexpr bool same(const T &p, const T &q) {
        if (!contains(p) || !contains(q)) return false;
        return get(p) == get(q);
    }

    // update [l, r] with value val
    // del: reflect effects of interval-delete
    // add: reflect effects of interval add
    template<class ADDFUNC, class DELFUNC> void update(T l, T r, const VAL &val, const ADDFUNC &add, const DELFUNC &del) {
        auto it = S.lower_bound(Node(l, 0, val));
        while (it != S.end() && it->l <= r) {
            if (it->l == r) {
                if (it->val ==val) {
                    del(r, it->r, val);
                    r = it->r;
                    it = S.erase(it);
                }
                break;
            }
            if (it->r <= r) {
                del(it->l, it->r, it->val);
                it = S.erase(it);
            } else {
                if (it->val == val) {
                    r = it->r;
                    del(it->l, it->r, it->val);
                    it = S.erase(it);
                } else {
                    del(it->l, r, it->val);
                    Node node = *it;
                    it = S.erase(it);
                    it = S.emplace_hint(it, r, node.r, node.val);
                }
            }
        }
        if (it != S.begin()) {
            it = prev(it);
            if (it->r == l) {
                if (it->val == val) {
                    del(it->l, it->r, it->val);
                    l = it->l;
                    it = S.erase(it);
                }
            } else if (l < it->r) {
                if (it->val == val) {
                    del(it->l, it->r, it->val);
                    l = min(l, it->l);
                    r = max(r, it->r);
                    it = S.erase(it);
                } else {
                    if (r < it->r) {
                        it = S.emplace_hint(next(it), r, it->r, it->val);
                        it = prev(it);
                    }
                    del(l, min(r, it->r), it->val);
                    Node node = *it;
                    it = S.erase(it);
                    it = S.emplace_hint(it, node.l, l, node.val);
                }
            }
        }
        if (it != S.end()) it = next(it);
        add(l, r, val);
        S.emplace_hint(it, l, r, val);
    }
    void update(const T &l, const T &r, const VAL &val) {
        update(l, r, val, [](T, T, VAL){}, [](T, T, VAL){});
    }
    template<class ADDFUNC, class DELFUNC> void update(T l, T r, const ADDFUNC &add, const DELFUNC &del) {
        update(l, r, VAL(), add, del);
    }
    void update(const T &l, const T &r) {
        update(l, r, VAL(), [](T, T, VAL){}, [](T, T, VAL){});
    }

    // remove
    template<class ADDFUNC, class DELFUNC> void remove(T l, T r, const ADDFUNC &add, const DELFUNC &del) {
        auto it = S.lower_bound(Node(l, 0, VAL()));
        while (it != S.end() && it->l <= r) {
            if (it->r <= r) {
                del(it->l, it->r, it->val);
                it = S.erase(it);
            } else {
                del(it->l, r, it->val);
                Node node = *it;
                it = S.erase(it);
                it = S.emplace_hint(it, r, node.r, node.val);
            }
        }
        if (it != S.begin()) {
            it = prev(it);
            if (l < it->r) {
                if (r < it->r) {
                    it = S.emplace_hint(next(it), r, it->r, it->val);
                    it = prev(it);
                }
                del(l, min(r, it->r), it->val);
                Node node = *it;
                it = S.erase(it);
                it = S.emplace_hint(it, node.l, l, node.val);
            }
        }
    }
    void remove(const T &l, const T &r) {
        remove(l, r, [](T, T, VAL){}, [](T, T, VAL){});
    }

    // debug
    friend ostream& operator << (ostream &s, const IntervalSet &ins) {
        for (auto e : ins.S) {
            s << "([" << e.l << ", " << e.r << "): " << e.val << ") ";
        }
        return s;
    }
};


//------------------------------//
// Examples
//------------------------------//

// 第六回 アルゴリズム実技検定 M - 等しい数
void PAST_6_M() {
    long long N, Q;
    cin >> N;
    vector<long long> A(N);
    map<long long, long long> cnt;
    for (int i = 0; i < N; i++) {
        cin >> A[i];
        cnt[A[i]]++;
    }
    long long res = 0;
    for (auto [val, num] : cnt) res += num * (num - 1) / 2;

    auto add = [&](int l, int r, long long val) -> void {
        long long before = cnt[val] * (cnt[val] - 1) / 2;
        cnt[val] += r - l;
        long long after = cnt[val] * (cnt[val] - 1) / 2;
        res += after - before;
    };
    auto del = [&](int l, int r, long long val) -> void {
        long long before = cnt[val] * (cnt[val] - 1) / 2;
        cnt[val] -= r - l;
        long long after = cnt[val] * (cnt[val] - 1) / 2;
        res += after - before;
    };
    IntervalSet<int, long long> ins(A);

    cin >> Q;
    while (Q--) {
        int l, r;
        long long val;
        cin >> l >> r >> val;
        l--;
        ins.update(l, r, val, add, del);
        cout << res << '\n';
    }
}

// RUPC 2018 G - Elevator
void RUPC_2018_G() {
    using fll = array<long long, 4>;
    int N, M, Q, t, l, r;
    cin >> N >> M >> Q;
    vector<fll> qs(M + Q);
    for (int i = 0; i < M; i++) {
        cin >> t >> l >> r, l--, r--;
        qs[i] = fll({-1, t*2+1, l*2, r*2+1});
    }
    for (int i = 0; i < Q; i++) {
        cin >> t >> l >> r, l--, r--;
        qs[i+M] = fll({i, t*2, l*2, r*2});
    }
    sort(qs.begin(), qs.end(), [&](const fll &p, const fll &q) { return p[1] < q[1]; });

    IntervalSet<long long> ins;
    vector<bool> res(Q, false);
    for (auto q : qs) {
        if (q[0] == -1) {
            ins.update(q[2], q[3]);
        } else {
            if (q[2] >= q[3] || ins.same(q[2], q[3])) res[q[0]] = true;
            else res[q[0]] = false;
        }
    }
    for (int q = 0; q < Q; ++q) cout << (res[q] ? "Yes" : "No") << endl;
}


int main() {
    //PAST_6_M();
    RUPC_2018_G();
}