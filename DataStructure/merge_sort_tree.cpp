//
// セグメント木を用いた merge-sort 木 (別名：領域木)
//
// verified:
//   ABC 339 G - Smaller Sum (merge-sort tree)
//     https://atcoder.jp/contests/abc339/tasks/abc339_g
//


#include <bits/stdc++.h>
using namespace std;


// Segment Tree
template<class Monoid> struct SegmentTree {
    using Func = function<Monoid(Monoid, Monoid)>;

    // core member
    int N;
    Func OP;
    Monoid IDENTITY;
    
    // inner data
    int log, offset;
    vector<Monoid> dat;

    // constructor
    SegmentTree() {}
    SegmentTree(int n, const Func &op, const Monoid &identity) {
        init(n, op, identity);
    }
    SegmentTree(const vector<Monoid> &v, const Func &op, const Monoid &identity) {
        init(v, op, identity);
    }
    void init(int n, const Func &op, const Monoid &identity) {
        N = n;
        OP = op;
        IDENTITY = identity;
        log = 0, offset = 1;
        while (offset < N) ++log, offset <<= 1;
        dat.assign(offset * 2, IDENTITY);
    }
    void init(const vector<Monoid> &v, const Func &op, const Monoid &identity) {
        init((int)v.size(), op, identity);
        build(v);
    }
    void pull(int k) {
        dat[k] = OP(dat[k * 2], dat[k * 2 + 1]);
    }
    void build(const vector<Monoid> &v) {
        assert(N == (int)v.size());
        for (int i = 0; i < N; ++i) dat[i + offset] = v[i];
        for (int k = offset - 1; k > 0; --k) pull(k);
    }
    int size() const {
        return N;
    }
    Monoid operator [] (int i) const {
        return dat[i + offset];
    }
    
    // update A[i], i is 0-indexed, O(log N)
    void set(int i, const Monoid &v) {
        assert(0 <= i && i < N);
        int k = i + offset;
        dat[k] = v;
        while (k >>= 1) pull(k);
    }
    
    // get [l, r), l and r are 0-indexed, O(log N)
    Monoid prod(int l, int r) {
        assert(0 <= l && l <= r && r <= N);
        Monoid val_left = IDENTITY, val_right = IDENTITY;
        l += offset, r += offset;
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) val_left = OP(val_left, dat[l++]);
            if (r & 1) val_right = OP(dat[--r], val_right);
        }
        return OP(val_left, val_right);
    }
    Monoid all_prod() {
        return dat[1];
    }
    
    // get max r such that f(v) = True (v = prod(l, r)), O(log N)
    // f(IDENTITY) need to be True
    int max_right(const function<bool(Monoid)> f, int l = 0) {
        if (l == N) return N;
        l += offset;
        Monoid sum = IDENTITY;
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(OP(sum, dat[l]))) {
                while (l < offset) {
                    l = l * 2;
                    if (f(OP(sum, dat[l]))) {
                        sum = OP(sum, dat[l]);
                        ++l;
                    }
                }
                return l - offset;
            }
            sum = OP(sum, dat[l]);
            ++l;
        } while ((l & -l) != l);  // stop if l = 2^e
        return N;
    }

    // get min l that f(get(l, r)) = True (0-indexed), O(log N)
    // f(IDENTITY) need to be True
    int min_left(const function<bool(Monoid)> f, int r = -1) {
        if (r == 0) return 0;
        if (r == -1) r = N;
        r += offset;
        Monoid sum = IDENTITY;
        do {
            --r;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(OP(dat[r], sum))) {
                while (r < offset) {
                    r = r * 2 + 1;
                    if (f(OP(dat[r], sum))) {
                        sum = OP(dat[r], sum);
                        --r;
                    }
                }
                return r + 1 - offset;
            }
            sum = OP(dat[r], sum);
        } while ((r & -r) != r);
        return 0;
    }
    
    // debug
    friend ostream& operator << (ostream &s, const SegmentTree &seg) {
        for (int i = 0; i < (int)seg.size(); ++i) {
            s << seg[i];
            if (i != (int)seg.size() - 1) s << " ";
        }
        return s;
    }
};


// ABC 339 G
void ABC_339_G() {
    int N;
    cin >> N;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    
    // merge sort tree (with 累積和)
    using Node = pair<vector<long long>, vector<long long>>;
    auto merge = [&](const vector<long long> &l, const vector<long long> &r) {
        vector<long long> res;
        int lp = 0, rp = 0;
        while (lp < l.size() || rp < r.size()) {
            if (lp == l.size()) res.push_back(r[rp++]);
            else if (rp == r.size()) res.push_back(l[lp++]);
            else if (l[lp] < r[rp]) res.push_back(l[lp++]);
            else res.push_back(r[rp++]);
        }
        return res;
    };
    auto merge2 = [&](const Node &l, const Node &r) {
        const auto &res = merge(l.first, r.first);
        vector<long long> sum({0});
        for (auto v : res) sum.push_back(sum.back() + v);
        return Node(res, sum);
    };
    Node unity(vector<long long>(), vector<long long>({0}));
    vector<Node> vec(N);
    for (int i = 0; i < N; ++i) {
        vec[i] = Node(vector<long long>({A[i]}), vector<long long>({0, A[i]}));
    }
    SegmentTree<Node> seg(vec, merge2, unity);
    
    // セグ木の構造を利用して X 以下の値の総和を求める
    auto getsum = [&](long long X, const Node &node) {
        int id = lower_bound(node.first.begin(), node.first.end(), X + 1) - node.first.begin();
        return node.second[id];
    };
    auto query = [&](auto query, int l, int r, long long X, int id, int left, int right) {
        if (l <= left && right <= r) {
            return getsum(X, seg.dat[id]);
        } else if (l < right && left < r) {
            int mid = (left + right) / 2;
            return query(query, l, r, X, id*2, left, mid)
                + query(query, l, r, X, id*2+1, mid, right);
        } else {
            return 0LL;
        }
    };
    
    // query
    int Q;
    cin >> Q;
    long long res = 0;
    while (Q--) {
        long long a, b, c;
        cin >> a >> b >> c;
        long long L = a ^ res, R = b ^ res, X = c ^ res;
        --L;
        res = query(query, L, R, X, 1, 0, seg.offset);
        cout << res << endl;
    }
}


int main() {
    ABC_339_G();
}
