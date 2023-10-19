//
// Wavelet Matrix (Preprocess: O(N log σ), Query: O(log σ))
//
// verified:
//   Yosupu Judge - Range Kth Smallest (for k_th_smallest)
//     https://judge.yosupo.jp/problem/range_kth_smallest
//
//   Yosupu Judge - Static Range Frequency (for rank)
//     https://judge.yosupo.jp/problem/static_range_frequency
//
//   Yosupu Judge - Rectangle Sum (for BIT on WM)
//     https://judge.yosupo.jp/problem/rectangle_sum
//
//   Yosupu Judge - Point Add Rectangle Sum (for BIT on WM)
//     https://judge.yosupo.jp/problem/point_add_rectangle_sum
//
//   AtCoder ABC 234 D - Prefix K-th Max (for k_th_largest)
//     https://atcoder.jp/contests/abc234/tasks/abc234_d
//
//   AtCoder ABC 281 E - Least Elements (for top_k_sum)
//     https://atcoder.jp/contests/abc281/tasks/abc281_e
//
//   AtCoder ABC 324 E - G - Generate Arrays (for range_freq)
//     https://atcoder.jp/contests/abc324/tasks/abc324_g
//
//   AOJ 1549 Hard Beans (for prev_value, next_value)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1549
//
//   AOJ 2426 Treasure Hunt (for BIT on WM)
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2426
//


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

// Wavelet Matrix (must vec[i] >= 0)
template<class T> struct WaveletMatrix {
    // inner data
    unsigned int n, height;
    vector<T> v;
    vector<BitVector> bv;
    vector<vector<long long>> sum;

    // constructor (sigma: the number of characters)
    WaveletMatrix() : n(0) {}
    WaveletMatrix(unsigned int n) : n(n), v(n) {}
    WaveletMatrix(const vector<T> &vec) : n(vec.size()), v(vec) {
        build();
    }
    void add(const T &val) {
        assert(v >= 0);
        v.push_back(v);
        n = v.size();
    }
    void set(unsigned int i, const T &val) {
        assert(i >= 0 && i < n && val >= 0);
        v[i] = val;
    }
    void build() {
        assert(n == (int)v.size());
        T mv = 1;
        for (int i = 0; i < n; ++i) mv = max(mv, v[i]);
        for (height = 1; mv != 0; mv >>= 1) ++height;
        vector<int> left(n), right(n), ord(n);
        iota(ord.begin(), ord.end(), 0);
        bv.assign(height, BitVector(n));
        sum.assign(height + 1, vector<long long>(n + 1, 0));
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
            for (int i = 0; i < n; ++i) sum[h][i + 1] = sum[h][i] + v[ord[i]];
        }
    }
    
    // access v[k]
    T access(int i) {
        T res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int i0 = bv[h].rank0(i);
            if (bv[h].get(i)) {
                i += bv[h].rank0() - i0;
                res |= T(1) << h;
            } else {
                i = i0;
            }
        }
        return res;
    }
    T operator [] (int i) {
        return access(i);
    }
    
    // count "i" s.t. v[i] = val, i \in [l, r)
    int rank(int l, int r, const T &val) {
        assert(0 <= l && l <= r && r <= n);
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if ((val >> h) & 1) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
            } else {
                l = l0;
                r = r0;
            }
        }
        return r - l;
    }
    
    // count "i" s.t. v[i] \in [lower, upper), i \in [l, r)
    int range_freq(int l, int r, const T &upper) {
        assert(0 <= l && l <= r && r <= n);
        int res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if ((upper >> h) & 1) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                res += r0 - l0;
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
    int range_freq(int l, int r, const T &lower, const T &upper) {
        return range_freq(l, r, upper) - range_freq(l, r, lower);
    }
    
    // the k-th (0-indexed) smallest value in [l, r)
    T k_th_smallest(int l, int r, int k) {
        assert(0 <= l && l <= r && r <= n);
        T res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if (r0 - l0 <= k) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                k -= r0 - l0;
                res |= T(1) << h;
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
    
    // the k-th (0-indexed) largest value in [l, r)
    T k_th_largest(int l, int r, int k) {
        assert(0 <= l && l <= r && r <= n);
       return k_th_smallest(l, r, r - l - k - 1);
    }
    
    // the sum of the top-k sum in [l, r)
    T top_k_sum(int l, int r, int k) {
        assert(0 <= l && l <= r && r <= n);
        if (l == r) return 0;
        T res = 0, val = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if (r0 - l0 <= k) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                k -= r0 - l0;
                val |= T(1) << h;
                res += sum[h][r0] - sum[h][l0];
            } else {
                l = l0;
                r = r0;
            }
        }
        res += val * k;
        return res;
    }
    
    // the max value (< val) in [l, r)
    T prev_value(int l, int r, T val) {
        assert(0 <= l && l <= r && r <= n);
        int num = range_freq(l, r, 0, val);
        if (num == 0) return T(-1);
        else return k_th_smallest(l, r, num - 1);
    }
    
    // the min value (>= val) in [l, r)
    T next_value(int l, int r, T val) {
        assert(0 <= l && l <= r && r <= n);
        int num = range_freq(l, r, 0, val);
        if (num == r - l) return T(-1);
        else return k_th_smallest(l, r, num);
    }
};


// 2D queries
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
    POS mi = 0;
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
        mi = min(mi, y);
    }
    int xid(POS x) const {
        return lower_bound(ps.begin(), ps.end(), Point(x, mi)) - ps.begin();
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
        bv.resize(height, BitVector(n));
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

void YosupoRangeKthSmallest() {
    int N, Q;
    cin >> N >> Q;
    vector<int> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    
    WaveletMatrix<int> wm(a);
    for (int q = 0; q < Q; ++q) {
        int l, r, k;
        cin >> l >> r >> k;
        cout << wm.k_th_smallest(l, r, k) << endl;
    }
}

void YosupoStaticRangeFrequency() {
    int N, Q;
    cin >> N >> Q;
    vector<int> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i];
    
    WaveletMatrix<int> wm(a);
    for (int q = 0; q < Q; ++q) {
        int l, r, x;
        cin >> l >> r >> x;
        cout << wm.rank(l, r, x) << endl;
    }
}

void YosupoRectangleSum() {
    BITonWaveletMatrix<int, long long> wm;
    
    int N, Q;
    cin >> N >> Q;
    vector<int> ix(N), iy(N), iw(N), lx(Q), ly(Q), rx(Q), ry(Q);
    for (int i = 0; i < N; ++i) {
        cin >> ix[i] >> iy[i] >> iw[i];
        wm.add_point(ix[i], iy[i]);
    }
    
    wm.build();
    
    for (int i = 0; i < N; ++i) {
        wm.add(ix[i], iy[i], iw[i]);
    }
    for (int q = 0; q < Q; ++q) {
        cin >> lx[q] >> ly[q] >> rx[q] >> ry[q];
        cout << wm.sum(lx[q], rx[q], ly[q], ry[q]) << endl;
    }
}

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

void ABC_234_D() {
    int N, K;
    cin >> N >> K;
    vector<int> P(N);
    for (int i = 0; i < N; ++i) cin >> P[i];
    
    WaveletMatrix wm(P);
    for (int i = K; i <= N; ++i) {
        cout << wm.k_th_largest(0, i, K - 1) << endl;
    }
}

void ABC_281_E() {
    int N, M, K;
    cin >> N >> M >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    
    WaveletMatrix wm(A);
    for (int i = 0; i < N - M + 1; ++i) {
        cout << wm.top_k_sum(i, i + M, K) << " ";
    }
    cout << endl;
}

void ABC_324_G() {
    int N, Q;
    cin >> N;
    vector<int> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    
    WaveletMatrix<int> wm(A);

    cin >> Q;
    vector<int> left(Q+1, 0), right(Q+1, N), lower(Q+1, 1), upper(Q+1, N+1);
    for (int q = 0; q < Q; ++q) {
        int t, s, v;
        cin >> t >> s >> v;
        if (t == 1) {
            // [left, x) 内部の lower 以上 upper 未満の値が v 個以上である最小の x を求める
            int low = left[s] - 1, high = N + 1;
            while (high - low > 1) {
                int mid = (low + high) / 2;
                if (wm.range_freq(left[s], mid, lower[s], upper[s]) >= v) high = mid;
                else low = mid;
            }
            
            // [left, high) と [high, right) とに分割する
            high = min(high, right[s]);
            left[q+1] = high;
            right[q+1] = right[s];
            lower[q+1] = lower[s];
            upper[q+1] = upper[s];
            right[s] = high;
        } else {
            // 値が v 未満と v 以上に分ける
            ++v;
            v = max(v, lower[s]);
            v = min(v, upper[s]);
            left[q+1] = left[s];
            right[q+1] = right[s];
            lower[q+1] = v;
            upper[q+1] = upper[s];
            upper[s] = v;
        }
        cout << wm.range_freq(left[q+1], right[q+1], lower[q+1], upper[q+1]) << endl;
    }
}

void AOJ_1549() {
    int N, Q;
    cin >> N;
    const int GETA = 1100000;
    vector<int> a(N);
    for (int i = 0; i < N; ++i) {
        cin >> a[i];
        a[i] += GETA;
    }
    
    WaveletMatrix<int> wm(a);
    cin >> Q;
    while (Q--) {
        int l, r, d;
        cin >> l >> r >> d;
        ++r;
        d += GETA;
        
        int nex = wm.next_value(l, r, d);
        int pre = wm.prev_value(l, r, d);
        int res = GETA*2;
        if (nex != -1) res = min(res, nex - d);
        if (pre != -1) res = min(res, d - pre);
        cout << res << endl;
    }
}

void AOJ_2426() {
    BITonWaveletMatrix<long long, long long> wm;
    
    int N, M;
    cin >> N >> M;
    vector<long long> x(N), y(N), lx(M), ly(M), rx(M), ry(M);
    for (int i = 0; i < N; ++i) {
        cin >> x[i] >> y[i];
        wm.add_point(x[i], y[i]);
    }
    
    wm.build();
    
    for (int i = 0; i < N; ++i) {
        wm.add(x[i], y[i], 1);
    }
    for (int q = 0; q < M; ++q) {
        cin >> lx[q] >> ly[q] >> rx[q] >> ry[q];
        ++rx[q], ++ry[q];
        cout << wm.sum(lx[q], rx[q], ly[q], ry[q]) << endl;
    }
}


int main() {
    //YosupoRangeKthSmallest();
    //YosupoStaticRangeFrequency();
    //YosupoRectangleSum();
    YosupoJudgePointAddRectangleSum();
    //ABC_234_D();
    //ABC_281_E();
    //ABC_324_G();
    //AOJ_1549();
    //AOJ_2426();
}

