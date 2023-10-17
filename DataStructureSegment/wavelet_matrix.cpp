//
// Wavelet Matrix (Preprocess: O(N log N), Query: O(log N))
//
// verified:
//   Yosupu Judge - Range Kth Smallest (for k_th_smallest)
//     https://judge.yosupo.jp/problem/range_kth_smallest
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
//   AOJ 2426 Treasure Hunt (for BIT on WM) 未 verify
//     https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=2426
//


#include <bits/stdc++.h>
using namespace std;



#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class T1, class T2> ostream& operator << (ostream &s, pair<T1,T2> P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, vector<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, deque<T> P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, vector<vector<T> > P)
{ for (int i = 0; i < P.size(); ++i) { s << endl << P[i]; } return s << endl; }
template<class T> ostream& operator << (ostream &s, set<T> P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T> ostream& operator << (ostream &s, multiset<T> P)
{ for (auto it : P) { s << "<" << it << "> "; } return s; }
template<class T1, class T2> ostream& operator << (ostream &s, map<T1,T2> P)
{ for (auto it : P) { s << "<" << it.first << "->" << it.second << "> "; } return s; }





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
    int height;
    vector<BitVector> bv;
    vector<vector<long long>> sum;

    // constructor (sigma: the number of characters)
    WaveletMatrix() {}
    WaveletMatrix(vector<T> vec) :
        WaveletMatrix(vec, *max_element(vec.begin(), vec.end()) + 1) {}
    WaveletMatrix(vector<T> vec, T sigma) {
        build(vec, sigma);
    }
    void build(vector<T> &vec, T sigma) {
        height = (sigma == 1) ? 1 : (64 - __builtin_clzll(sigma - 1));
        bv.resize(height);
        vector<T> A = vec;
        sum.resize(height + 1);
        for (int h = 0; h < height; ++h) {
            bv[h].resize(vec.size());
            for (int j = 0; j < vec.size(); ++j) {
                bv[h].set(j, get(vec[j], height - h - 1));
            }
            bv[h].build();
            stable_partition(vec.begin(), vec.end(), [&](int c) {
                return !get(c, height - h - 1);
            });
        }
        for (int h = 0; h <= height; ++h) {
            sum[h].resize((int)A.size() + 1);
            for (int j = 1; j <= (int)A.size(); ++j) {
                sum[h][j] = sum[h][j - 1] + A[j - 1];
            }
            if (h == height) break;
            stable_partition(A.begin(), A.end(), [&](int c) {
                return !get(c, height - h - 1);
            });
        }
    }
    
    // the i-th bit of "val" (0 or 1)
    int get(const T val, const int i) {
        return val >> i & 1;
    }
    
    // the number of "val" in [l, r)
    int rank(const T val, const int l, const int r) {
        return rank(val, r) - rank(val, l);
    }
    
    // the number of "val" in [0, i)
    int rank(T val, int i) {
        int p = 0;
        for (int h = 0; h < height; ++h) {
            int p0 = bv[h].rank0(p), i0 = bv[h].rank0(i);
            if (get(val, height - h - 1)) {
                p += bv[h].rank0() - p0;
                i += bv[h].rank0() - i0;
            } else {
                p = p0;
                i = i0;
            }
        }
        return i - p;
    }
    
    // the k-th (0-indexed) smallest value in [l, r)
    T k_th_smallest(int k, int l, int r) {
        T res = 0;
        for (int h = 0; h < height; ++h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if (r0 - l0 > k) {
                l = l0;
                r = r0;
            } else {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                k -= r0 - l0;
                res |= (1LL << (height - h - 1));
            }
        }
        return res;
    }
    
    // the k-th (0-indexed) largest value in [l, r)
    T k_th_largest(int k, int l, int r) {
       return k_th_smallest(r - l - k - 1, l, r);
    }
    
    // the sum of the top-k sum in [l, r)
    T top_k_sum(int k, int l, int r) {
        if (l == r) return 0;
        T res = 0, val = 0;
        for (int h = 0; h < height; ++h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if (r0 - l0 > k) {
                l = l0;
                r = r0;
            } else {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                k -= r0 - l0;
                val |= (1LL << (height - h - 1));
                res += sum[h + 1][r0] - sum[h + 1][l0];
            }
        }
        res += (long long)val * k;
        return res;
    }
    
    // the number of value between [loewr, upper) in [l, r)
    /*
    int range_freq(int l, int r, T upper) {
        if (upper >= (1LL << height)) return r - l;
        int res = 0;
        for (int h = height - 1; h >= 0; --h) {
            int l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if ((upper >> h) & 1) {
                l += bv[h].rank0() - l0;
                r += bv[h].rank0() - r0;
                res += r - l;
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
    int range_freq(int l, int r, T lower, T upper) {
        return range_freq(l, r, upper) - range_freq(l, r, lower);
    }
     */
    int range_freq(const int i, const int j, const T lower, const T upper,
                   const int l, const int r, const int h) {
        if (i == j || r <= lower || upper <= l) return 0;
        int mid = (l + r) >> 1;
        if (lower <= l && r <= upper) {
            return j - i;
        } else {
            int i0 = bv[h].rank0(i), j0 = bv[h].rank0(j);
            T left = range_freq(i0, j0, lower, upper, l, mid, h + 1);
            T right = range_freq(i + bv[h].rank0() - i0, j + bv[h].rank0() - j0,
                                 lower, upper, mid, r, h + 1);
            return left + right;
        }
    }
    int range_freq(const int l, const int r, const T lower, const T upper) {
        return range_freq(l, r, lower, upper, 0, 1LL << height, 0);
    }
    
    // the minmum value between [lower, upper) in [l, r)
    T range_min(const int i, const int j, const T lower, const T upper,
                  const int l, const int r, const int h, const T val) {
        if (i == j || r <= lower || upper <= l) return -1;
        if (r - l == 1) return val;
        int mid = (l + r) >> 1;
        int i0 = bv[h].rank0(i), j0 = bv[h].rank0(j);
        T res = range_min(i0, j0, lower, upper, l, mid, h + 1, val);
        if (res < 0) {
            return range_min(i + bv[h].rank0() - i0, j + bv[h].rank0() - j0,
                             lower, upper, mid, r, h + 1, val + (1LL << (height - h - 1)));
        } else {
            return res;
        }
    }
    T range_min(int l, int r, T lower, T upper) {
        return range_min(l, r, lower, upper, 0, 1LL << height, 0, 0);
    }
    
    // the max value (< val) in [l, r)
    T prev_value(int l, int r, T val) {
        int num = range_freq(l, r, 0, val);
        if (num == 0) return T(-1);
        else return k_th_smallest(num - 1, l, r);
    }
    
    // the min value (>= val) in [l, r)
    T next_value(int l, int r, T val) {
        int num = range_freq(l, r, 0, val);
        if (num == r - l) return T(-1);
        else return k_th_smallest(num, l, r);
    }
};

// 2D range count
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
    // add_point() cannot be used after init()
    BITonWaveletMatrix() {}
    BITonWaveletMatrix(vector<Point> vec) {
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
        for (height = 0; mv != 0; mv >>= 1) ++height;
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
        cout << wm.k_th_smallest(k, l, r) << endl;
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
        cout << wm.k_th_largest(K - 1, 0, i) << endl;
    }
}

void ABC_281_E() {
    int N, M, K;
    cin >> N >> M >> K;
    vector<long long> A(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    
    WaveletMatrix wm(A);
    for (int i = 0; i < N - M + 1; ++i) {
        cout << wm.top_k_sum(K, i, i + M) << " ";
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
    for (int i = 0; i < M; ++i) {
        cin >> lx[i] >> ly[i] >> rx[i] >> ry[i];
        ++rx[i], ++ry[i];
    }
    
    wm.build();
    
    for (int i = 0; i < N; ++i) {
        wm.add(x[i], y[i], 1);
    }
    for (int i = 0; i < M; ++i) {
        cout << wm.sum(lx[i], rx[i], ly[i], ry[i]) << endl;
    }
}


int main() {
    //YosupoRangeKthSmallest();
    YosupoJudgePointAddRectangleSum();
    //ABC_234_D();
    //ABC_281_E();
    //ABC_324_G();
    //AOJ_1549();
    //AOJ_2426();
}

