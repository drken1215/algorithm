//
// BIT on Wavelet Matrix (Preprocess: O(N log N), Query: O((log N)^2))
//
// verified:
//   Yosupu Judge - Point Add Rectangle Sum (for BIT on WM)
//     https://judge.yosupo.jp/problem/point_add_rectangle_sum
//

/*
 　クエリ先読みを前提として、一点加算のみ可能にしたウェーブレット行列
 　　・加算クエリが発生する座標 (x, y) をあらかじめ wm.add_point(x, y) を用いて登録する
 　　・取得クエリが発生する区間 (lx, rx, ly, ry) の登録は不要 (内部で自動的に座標圧縮される)
 　　・登録後に wm.build() する (以後、add_point(x, y) は使用不可)
 
 　　・その後は、以下のクエリを O(log N) で実行
 　　　　・一点加算クエリ add(x, y, w)
 　　　　・区間総和取得クエリ sum(lx, rx, ly, ry)
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


// BIT on Wavelet Matrix
template<class POS, class VAL> struct BITonWaveletMatrix {
    // inner data
    struct BIT {
        VAL UNITY_SUM = 0;
        int N;
        vector<VAL> dat;
        
        // [0, n)
        BIT() {}
        BIT(int n, VAL unity = 0) : UNITY_SUM(unity), N(n), dat(n, unity) { }
        void init(int n) {
            N = n;
            dat.assign(n, UNITY_SUM);
        }
        
        // a is 0-indexed
        void add(int a, VAL x) {
            for (int i = a; i < (int)dat.size(); i |= i + 1)
                dat[i] = dat[i] + x;
        }
        
        // [0, a), a is 0-indexed
        VAL sum(int a) {
            VAL res = UNITY_SUM;
            for (int i = a - 1; i >= 0; i = (i & (i + 1)) - 1)
                res = res + dat[i];
            return res;
        }
        
        // [a, b), a and b are 0-indexed
        VAL sum(int a, int b) {
            return sum(b) - sum(a);
        }
        
        // debug
        void print() {
            for (int i = 0; i < (int)dat.size(); ++i)
                cout << sum(i, i + 1) << ",";
            cout << endl;
        }
    };
    
    using Point = pair<POS, POS>;
    int n, height;
    vector<BitVector> bv;
    vector<Point> ps;
    vector<POS> ys;
    vector<BIT> bit;

    // constructor (sigma: the number of characters)
    // add_point() cannot be used after build()
    BITonWaveletMatrix() {}
    BITonWaveletMatrix(const vector<Point> &vec) {
        for (auto [x, y] : vec) add_point(x, y);
        build();
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
    void build() {
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
        bit.assign(height + 1, BIT(n));
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
    
    // add
    void add(const POS x, const POS y, const VAL val) {
        int i = lower_bound(ps.begin(), ps.end(), Point(x, y)) - ps.begin();
        int j = yid(y);
        for (int h = height - 1; h >= 0; --h) {
            int i0 = bv[h].rank0(i);
            if ((j >> h) & 1) {
                i += bv[h].rank0() - i0;
            } else {
                i = i0;
            }
            bit[h].add(i, val);
        }
    }
    
    // sum
    VAL inner_sum(int l, int r, const POS upper) {
        assert(0 <= l && r <= n);
        VAL res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if ((upper >> h) & 1) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                res += bit[h].sum(l0, r0);
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
    VAL sum(const POS lx, const POS rx, const POS ly, const POS ry) {
        int l = xid(lx), r = xid(rx);
        return inner_sum(l, r, yid(ry)) - inner_sum(l, r, yid(ly));
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

void YosupoJudgePointAddRectangleSum() {
    BITonWaveletMatrix<int, long long> wm;
    
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
    
    wm.build();
    for (int i = 0; i < N; ++i) {
        wm.add(ix[i], iy[i], iw[i]);
    }
    for (int q = 0; q < Q; ++q) {
        if (type[q] == 0) {
            wm.add(lx[q], ly[q], rx[q]);
        } else {
            cout << wm.sum(lx[q], rx[q], ly[q], ry[q]) << endl;
        }
    }
}


int main() {
    YosupoJudgePointAddRectangleSum();
}
