//
// Wavelet Matrix (Preprocess: O(N log N), Query: O(log N))
//
// verified:
//   AtCoder ABC 234 D - Prefix K-th Max (for k_th_largest)
//     https://atcoder.jp/contests/abc281/tasks/abc281_e
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


#include <bits/stdc++.h>
using namespace std;


// Bit Vector
struct BitRank {
    // block: bit vector
    // count: the number of 1 within each block
    vector<unsigned long long> block;
    vector<unsigned int> count;
    
    // constructor
    BitRank() {}
    void resize(const unsigned int num) {
        block.resize(((num + 1) >> 6) + 1, 0);
        count.resize(block.size(), 0);
    }
    
    // set val(0 or 1) onto i-th bit
    void set(const unsigned int i, const unsigned long long val) {
        block[i >> 6] |= (val << (i & 63));
    }
    void build() {
        for (unsigned int i = 1; i < block.size(); i++) {
            count[i] = count[i - 1] + __builtin_popcountll(block[i - 1]);
        }
    }
    
    // the number of 1 in [0, i)
    unsigned int rank1(const unsigned int i) const {
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
};

// Wavelet Matrix
template<class T> struct WaveletMatrix {
    // inner data
    int height;
    vector<BitRank> bv;
    vector<int> pos;
    vector<vector<long long>> rui;

    // constructor (sigma: the number of characters)
    WaveletMatrix() {}
    WaveletMatrix(vector<T> vec) :
        WaveletMatrix(vec, *max_element(vec.begin(), vec.end()) + 1) {}
    WaveletMatrix(vector<T> vec, T sigma) {
        init(vec, sigma);
    }
    void init(vector<T> &vec, T sigma) {
        height = (sigma == 1) ? 1 : (64 - __builtin_clzll(sigma - 1));
        bv.resize(height), pos.resize(height);
        vector<T> A = vec;
        rui.resize(height + 1);
        for (int i = 0; i < height; ++i) {
            bv[i].resize(vec.size());
            for (int j = 0; j < vec.size(); ++j) {
                bv[i].set(j, get(vec[j], height - i - 1));
            }
            bv[i].build();
            auto it = stable_partition(vec.begin(), vec.end(), [&](int c) {
                return !get(c, height - i - 1);
            });
            pos[i] = it - vec.begin();
        }
        for (int i = 0; i <= height; ++i) {
            rui[i].resize((int)A.size() + 1);
            for (int j = 1; j <= (int)A.size(); ++j) {
                rui[i][j] = rui[i][j - 1] + A[j - 1];
            }
            if (i == height) break;
            stable_partition(A.begin(), A.end(), [&](int c) {
                return !get(c, height - i - 1);
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
        for (int j = 0; j < height; ++j) {
            if (get(val, height - j - 1)) {
                p = pos[j] + bv[j].rank1(p);
                i = pos[j] + bv[j].rank1(i);
            } else {
                p = bv[j].rank0(p);
                i = bv[j].rank0(i);
            }
        }
        return i - p;
    }
    
    // the k-th (0-indexed) smallest value in [l, r)
    T k_th_smallest(int k, int l, int r) {
        T res = 0;
        for (int i = 0; i < height; ++i) {
            const int j = bv[i].rank0(l, r);
            if (j > k){
                l = bv[i].rank0(l);
                r = bv[i].rank0(r);
            } else {
                l = pos[i] + bv[i].rank1(l);
                r = pos[i] + bv[i].rank1(r);
                k -= j;
                res |= (1LL << (height - i - 1));
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
        for (int i = 0; i < height; ++i) {
            const int j = bv[i].rank0(l, r);
            if (j > k) {
                l = bv[i].rank0(l);
                r = bv[i].rank0(r);
            } else {
                int l2 = bv[i].rank0(l);
                int r2 = bv[i].rank0(r);
                res += rui[i + 1][r2] - rui[i + 1][l2];
                l = pos[i] + bv[i].rank1(l);
                r = pos[i] + bv[i].rank1(r);
                k -= j;
                val |= (1LL << (height - i - 1));
            }
        }
        res += (long long)val * k;
        return res;
    }
    
    // the number of value between [loewr, upper) in [l, r)
    int range_freq(const int i, const int j, const T lower, const T upper,
                   const int l, const int r, const int x) {
        if (i == j || r <= lower || upper <= l) return 0;
        int mid = (l + r) >> 1;
        if (lower <= l && r <= upper) {
            return j - i;
        } else {
            T left = range_freq(bv[x].rank0(i), bv[x].rank0(j), lower, upper, l, mid, x + 1);
            T right = range_freq(pos[x] + bv[x].rank1(i), pos[x] + bv[x].rank1(j),
                                 lower, upper, mid, r, x + 1);
            return left + right;
        }
    }
    int range_freq(const int l, const int r, const T lower, const T upper) {
        return range_freq(l, r, lower, upper, 0, 1LL << height, 0);
    }
    
    // the minmum value between [lower, upper) in [l, r)
    T range_min(const int i, const int j, const T lower, const T upper,
                  const int l, const int r, const int x, const T val) {
        if (i == j || r <= lower || upper <= l) return -1;
        if (r - l == 1) return val;
        int mid = (l + r) >> 1;
        T res = range_min(bv[x].rank0(i), bv[x].rank0(j), lower, upper, l, mid, x + 1, val);
        if (res < 0) {
            return range_min(pos[x] + bv[x].rank1(i), pos[x] + bv[x].rank1(j),
                             lower, upper, mid, r, x + 1, val + (1LL << (height - x - 1)));
        } else return res;
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
template<class T> struct OrthogonalRangeCount {
    // inner data
    using ptt = pair<T, T>;
    const int sz;
    vector<T> X, Y;
    WaveletMatrix<T> wm;
    
    // constructor
    OrthogonalRangeCount(vector<ptt> candidate) : sz((int)candidate.size()), X(sz), Y(sz) {
        sort(candidate.begin(), candidate.end());
        vector<int> vec(sz);
        for (int i = 0; i < sz; ++i) {
            X[i] = candidate[i].first, Y[i] = candidate[i].second;
        }
        sort(Y.begin(), Y.end());
        Y.erase(unique(Y.begin(), Y.end()), Y.end());
        for (int i = 0; i < sz; ++i) {
            vec[i] = lower_bound(Y.begin(), Y.end(), candidate[i].second) - Y.begin();
        }
        wm.init(vec, Y.size());
    }
    
    // the number of points in [lx, rx) × [ly, ry)
    int query(const T lx, const T ly, const T rx, const T ry){
        const int lxid = lower_bound(X.begin(), X.end(), lx) - X.begin();
        const int rxid = lower_bound(X.begin(), X.end(), rx) - X.begin();
        const int lyid = lower_bound(Y.begin(), Y.end(), ly) - Y.begin();
        const int ryid = lower_bound(Y.begin(), Y.end(), ry) - Y.begin();
        if (lxid >= rxid || lyid >= ryid) return 0;
        return wm.range_freq(lxid, rxid, lyid, ryid);
    }
};



/*/////////////////////////////*/
// Examples
/*/////////////////////////////*/

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


int main() {
    //ABC_234_D();
    //ABC_281_E();
    //ABC_324_G();
    AOJ_1549();
}

