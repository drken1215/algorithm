//
// セグメント木 on Wavelet Matrix (Preprocess: O(N log N), Query: O((log N)^2))
//
// verified:
//   Yosupu Judge - Point Add Rectangle Sum
//     https://judge.yosupo.jp/problem/point_add_rectangle_sum
//

/*
 　クエリ先読みを前提として、一点更新のみ可能にしたウェーブレット行列
 　　・更新クエリが発生する座標 (x, y) をあらかじめ wm.add_point(x, y) を用いて登録する
 　　・取得クエリが発生する区間 (lx, rx, ly, ry) の登録は不要 (内部で自動的に座標圧縮される)
 　　・登録後に wm.build() する (以後、add_point(x, y) は使用不可)
 
 　　・その後は、以下のクエリを O(log N) で実行
 　　　　・一点更新クエリ set(x, y, w)
 　　　　・区間取得クエリ prod(lx, rx, ly, ry)
 
 　　セグメント木に載せられるデータ構造は「アーベル群」、可換かつ逆元をもつ
 */


#include <bits/stdc++.h>
using namespace std;


// Bit Vector (for 64-bit non-negative integer)
struct BitVector {
    // block: bit vector
    // count: the number of 1 within each block
    unsigned int n, zeros;
    vector<unsigned long long> block;
    vector<unsigned int> count;
    
    // constructor
    BitVector() {}
    BitVector(const unsigned int num) {
        resize(num);
    }
    void resize(const unsigned int num) {
        n = num;
        block.assign(((num + 1) >> 6) + 1, 0);
        count.assign(block.size(), 0);
    }
    
    // set val(0 or 1) onto i-th bit, get i-th bit of val(0 or 1)
    void set(const unsigned int i, const unsigned long long val = 1LL) {
        assert((i >> 6) < block.size());
        block[i >> 6] |= (val << (i & 63));
    }
    unsigned int get(const unsigned int i) const {
        assert((i >> 6) < block.size());
        return (const unsigned int)(block[i >> 6] >> (i & 63)) & 1;
    }
    void build() {
        for (unsigned int i = 1; i < block.size(); i++) {
            count[i] = count[i - 1] + __builtin_popcountll(block[i - 1]);
        }
        zeros = rank0(n);
    }
    
    // the number of 1 in [0, i)
    unsigned int rank1(const unsigned int i) const {
        assert((i >> 6) < count.size());
        assert((i >> 6) < block.size());
        return count[i >> 6] +
        __builtin_popcountll(block[i >> 6] & ((1ULL << (i & 63)) - 1ULL));
    }
    // the number of 1 in [i, j)
    unsigned int rank1(const unsigned int i, const unsigned int j) const {
        return rank1(j) - rank1(i);
    }
    // the number of 0 in [0, i)
    unsigned int rank0(const unsigned int i) const {
        return i - rank1(i);
    }
    // the number of 0 in [i, j)
    unsigned int rank0(const unsigned int i, const unsigned int j) const {
        return rank0(j) - rank0(i);
    }
    // the number of 0 in [0, n)
    unsigned int rank0() const {
        return zeros;
    }
};

// Segment Tree on Wavelet Matrix
template<class POS, class Abel> struct SegmentTreeOnWaveletMatrix {
    using Func = function<Abel(Abel, Abel)>;
    using InvFunc = function<Abel(Abel)>;
    struct SegmentTree {
        // core member
        int N;
        Func OP;
        Abel IDENTITY;
        
        // inner data
        int log, offset;
        vector<Abel> dat;

        // constructor
        SegmentTree() {}
        SegmentTree(int n, const Func &op, const Abel &identity) {
            init(n, op, identity);
        }
        void init(int n, const Func &op, const Abel &identity) {
            N = n;
            OP = op;
            IDENTITY = identity;
            log = 0, offset = 1;
            while (offset < N) ++log, offset <<= 1;
            dat.assign(offset * 2, IDENTITY);
        }
        void pull(int k) {
            dat[k] = OP(dat[k * 2], dat[k * 2 + 1]);
        }
        
        // update A[i], i is 0-indexed, O(log N)
        void set(int i, const Abel &v) {
            assert(0 <= i && i < N);
            int k = i + offset;
            dat[k] = v;
            while (k >>= 1) pull(k);
        }
        
        // get [l, r), l and r are 0-indexed, O(log N)
        Abel prod(int l, int r) {
            assert(0 <= l && l <= r && r <= N);
            Abel val_left = IDENTITY, val_right = IDENTITY;
            l += offset, r += offset;
            for (; l < r; l >>= 1, r >>= 1) {
                if (l & 1) val_left = OP(val_left, dat[l++]);
                if (r & 1) val_right = OP(dat[--r], val_right);
            }
            return OP(val_left, val_right);
        }
    };
    
    using Point = pair<POS, POS>;
    Func OP;
    InvFunc IOP;
    Abel IDENTITY;
    int n, height;
    vector<BitVector> bv;
    vector<Point> ps;
    vector<POS> ys;
    vector<SegmentTree> seg;

    // constructor (sigma: the number of characters)
    // add_point() cannot be used after build()
    SegmentTreeOnWaveletMatrix() {}
    SegmentTreeOnWaveletMatrix(const vector<Point> &vec) {
        for (auto [x, y] : vec) add_point(x, y);
    }
    SegmentTreeOnWaveletMatrix(const vector<Point> &vec,
                               const Func &op, const InvFunc &iop, const Abel &identity) {
        for (auto [x, y] : vec) add_point(x, y);
        build(op, iop, identity);
    }
    void add_point(POS x, POS y) {
        ps.emplace_back(x, y);
        ys.emplace_back(y);
    }
    int xid(POS x) const {
        return lower_bound(ps.begin(), ps.end(), Point(x, 0)) - ps.begin();
    }
    int yid(POS y) const {
        return lower_bound(ys.begin(), ys.end(), y) - ys.begin();
    }
    void build(const Func &op, const InvFunc &iop, const Abel &identity) {
        OP = op, IOP = iop, IDENTITY = identity;
        sort(ps.begin(), ps.end());
        ps.erase(unique(ps.begin(), ps.end()), ps.end());
        n = (int)ps.size();
        sort(ys.begin(), ys.end());
        ys.erase(unique(ys.begin(), ys.end()), ys.end());
        vector<int> v(n), left(n), right(n), ord(n);
        int mv = 1;
        for (int i = 0; i < n; ++i) {
            v[i] = yid(ps[i].second);
            mv = max(mv, v[i]);
        }
        for (height = 1; mv != 0; mv >>= 1) ++height;
        iota(ord.begin(), ord.end(), 0);
        bv.assign(height, BitVector(n));
        seg.assign(height + 1, SegmentTree(n, op, identity));
        for (int h = height - 1; h >= 0; --h) {
            int l = 0, r = 0;
            for (int i = 0; i < n; ++i) {
                if ((v[ord[i]] >> h) & 1) {
                    bv[h].set(i);
                    right[r++] = ord[i];
                } else {
                    left[l++] = ord[i];
                }
            }
            bv[h].build();
            ord.swap(left);
            for (int i = 0; i < r; ++i) ord[i + l] = right[i];
        }
    }
    
    // set
    void set(const POS x, const POS y, const Abel val) {
        int i = lower_bound(ps.begin(), ps.end(), Point(x, y)) - ps.begin();
        int j = yid(y);
        for (int h = height - 1; h >= 0; --h) {
            int i0 = bv[h].rank0(i);
            if ((j >> h) & 1) {
                i += bv[h].rank0() - i0;
            } else {
                i = i0;
            }
            seg[h].set(i, val);
        }
    }
    
    // prod
    Abel inner_prod(int l, int r, int upper) {
        assert(0 <= l && r <= n);
        Abel res = IDENTITY;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if ((upper >> h) & 1) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                res = OP(res, seg[h].prod(l0, r0));
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
    Abel prod(const POS lx, const POS rx, const POS ly, const POS ry) {
        int l = xid(lx), r = xid(rx);
        return OP(inner_prod(l, r, yid(ry)), IOP(inner_prod(l, r, yid(ly))));
    }
};



//------------------------------//
// Examples
//------------------------------//

void YosupoJudgePointAddRectangleSum() {
    SegmentTreeOnWaveletMatrix<int, long long> wm;
    
    // monoid function of Segment Tree
    auto op = [&](long long l, long long r) { return l + r; };
    auto iop = [&](long long v) { return -v; };
    long long identity = 0;
    
    int N, Q;
    cin >> N >> Q;
    vector<int> ix(N), iy(N), iw(N), type(Q), lx(Q), ly(Q), rx(Q), ry(Q);
    for (int i = 0; i < N; ++i) {
        cin >> ix[i] >> iy[i] >> iw[i];
        wm.add_point(ix[i], iy[i]);
    }
    for (int q = 0; q < Q; ++q) {
        cin >> type[q];
        if (type[q] == 0) {
            cin >> lx[q] >> ly[q] >> rx[q];
            wm.add_point(lx[q], ly[q]);
        }
        else cin >> lx[q] >> ly[q] >> rx[q] >> ry[q];
    }
    
    wm.build(op, iop, identity);
    for (int i = 0; i < N; ++i) {
        long long original_val = wm.prod(ix[i], ix[i]+1, iy[i], iy[i]+1);
        wm.set(ix[i], iy[i], op(original_val, iw[i]));
    }
    for (int q = 0; q < Q; ++q) {
        if (type[q] == 0) {
            long long original_val = wm.prod(lx[q], lx[q]+1, ly[q], ly[q]+1);
            wm.set(lx[q], ly[q], op(original_val, rx[q]));
        } else {
            cout << wm.prod(lx[q], rx[q], ly[q], ry[q]) << endl;
        }
    }
}


int main() {
    YosupoJudgePointAddRectangleSum();
}

