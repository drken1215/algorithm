//
// binary 行列 (with bitset 高速化)
//
// verified:
//   Yosupo Library Checker - Matrix Product (Mod 2)
//     https://judge.yosupo.jp/problem/matrix_product_mod_2
//
//   Yosupo Library Checker - Determinant of Matrix (Mod 2)
//     https://judge.yosupo.jp/problem/matrix_det_mod_2
//
//   Yosupo Library Checker - Rank of Matrix (Mod 2)
//     https://judge.yosupo.jp/problem/matrix_rank_mod_2
//
//   Yosupo Library Checker - System of Linear Equations (Mod 2)
//     https://judge.yosupo.jp/problem/system_of_linear_equations_mod_2
//
//   Yosupo Library Checker - Inverse Matrix (Mod 2)
//     https://judge.yosupo.jp/problem/inverse_matrix_mod_2
//


#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Utility
//------------------------------//

// output
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template<class S, class T> ostream& operator << (ostream &s, const pair<S, T> &P)
{ return s << '<' << P.first << ", " << P.second << '>'; }
template<class T> ostream& operator << (ostream &s, const vector<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }
template<class T> ostream& operator << (ostream &s, const deque<T> &P)
{ for (int i = 0; i < P.size(); ++i) { if (i > 0) { s << " "; } s << P[i]; } return s; }

// num of i such that (x & (1 << i)) != 0
int popcnt(int x) { return __builtin_popcount(x); }
int popcnt(unsigned int x) { return __builtin_popcount(x); }
int popcnt(long long x) { return __builtin_popcountll(x); }
int popcnt(unsigned long long x) { return __builtin_popcountll(x); }
int popcnt_mod2(int x) { return __builtin_parity(x); }
int popcnt_mod2(unsigned int x) { return __builtin_parity(x); }
int popcnt_mod2(long long x) { return __builtin_parityll(x); }
int popcnt_mod2(unsigned long long x) { return __builtin_parityll(x); }

// 動的 bitset
struct DynamicBitset {
    using u64 = unsigned long long;
    constexpr int topbit(u64 x) const { return (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
    constexpr int lowbit(u64 x) const { return (x == 0 ? -1 : __builtin_ctzll(x)); }
    static string CONV[256];  // for to_string()
    
    // inner values
    int N;
    vector<u64> dat;

    // constructors
    DynamicBitset(int N = 0, int x = 0) : N(N) {
        assert(x == 0 || x == 1);
        u64 v = (x == 0 ? 0 : -1);
        dat.assign((N + 63) >> 6, v);
        if (N) dat.back() >>= (64 * dat.size() - N);
    }
    DynamicBitset(const DynamicBitset&) = default;
    DynamicBitset& operator = (const DynamicBitset&) = default;
    DynamicBitset(const string &S) : N((int)S.size()) {
        assign(N, 0);
        for (int i = 0; i < N; i++) {
            if (S[i] == '1') dat[i >> 6] |= u64(1) << (i & 63);
        }
    }
    DynamicBitset(const vector<int> &v) : N((int)v.size()) {
        assign(N, 0);
        for (int i = 0; i < N; i++) {
            if (v[i] == 1) dat[i >> 6] |= u64(1) << (i & 63);
        }
    }
    constexpr int size() const { return N; }
    void resize(int siz) {
        dat.resize((siz + 63) >> 6);
        int rem = siz & 63;
        if (rem != 0) {
            u64 mask = (u64(1) << rem) - 1;
            dat.back() &= mask;
        }
        N = siz;
    }
    void assign(int siz, int x = 0) {
        assert(x == 0 || x == 1);
        N = siz;
        u64 v = (x == 0 ? 0 : -1);
        dat.assign((N + 63) >> 6, v);
        if (N) dat.back() >>= (64 * dat.size() - N);
    }
    int back() const {
        assert(N > 0);
        return (*this)[N - 1];
    }
    void push_back(bool v) {
        resize(N + 1);
        (*this)[size() - 1] = v;
    }
    static void pre_to_string() {
        for (int s = 0; s < 256; s++) {
            string x;
            for (int i = 0; i < 8; i++) x += '0' + (s >> i & 1);
            CONV[s] = x;
        }
    }
    string to_string() const {
        if (CONV[0].empty()) pre_to_string();
        string S;
        for (auto &x : dat) { 
            for (int i = 0; i < 8; i++) S += CONV[(x >> (8 * i) & 255)];
        }
        S.resize(N);
        return S;
    }
    
    // getters
    struct Proxy {
        int index;
        vector<u64> &dat;
        Proxy(vector<u64> &d, int i) : dat(d), index(i) {}
        operator bool () const {
            return (dat[index >> 6] >> (index & 63)) & 1;
        }
        Proxy& operator = (const Proxy& v) {
            return (*this) = bool(v);
        }
        Proxy& operator = (u64 value) {
            dat[index >> 6] &= ~(u64(1) << (index & 63));
            dat[index >> 6] |= (value & 1) << (index & 63);
            return *this;
        }
        void flip() {
            dat[index >> 6] ^= (u64(1) << (index & 63));
        }
    };
    Proxy operator [] (int i) { return Proxy(dat, i); }
    bool operator [] (int i) const { return (dat[i >> 6] >> (i & 63)) & 1; }
    bool any() const {  // [0, N)
        for (auto val : dat) if (val) return true;
        return false;
    }
    int count() const {  // [0, N)
        int res = 0;
        for (auto val : dat) res += popcnt(val);
        return res;
    }
    int count(int L, int R) const {  // [L, R)
        assert(L <= R);
        int res = 0;
        while ((L < R) && (L & 63)) res += (*this)[L++];
        while ((L < R) && (R & 63)) res += (*this)[--R];
        for (int i = (L >> 6); i < (R >> 6); i++) res += popcnt(dat[i]);
        return res;
    }
    int dot(const DynamicBitset &r) const {
        assert(size() == r.size());
        int res = 0;
        for (int i = 0; i < (int)dat.size(); i++) res += popcnt(dat[i] & r.dat[i]);
        return res;
    }
    int dot_mod2(const DynamicBitset &r) const {
        assert(size() == r.size());
        int res = 0;
        for (int i = 0; i < (int)dat.size(); i++) res ^= popcnt_mod2(dat[i] & r.dat[i]);
        return res;
    }
    int next(int val) const {  // (include val)
        if (val < 0) val = 0;
        if (val >= N) return N;
        int k = val >> 6;
        {
            u64 x = dat[val >> 6];
            int s = val & 63;
            x = (x >> s) << s;
            if (x) return (k << 6) | lowbit(x);
        }
        for (int i = k + 1; i < (int)dat.size(); i++) {
            if (dat[i] == 0) continue;
            return (i << 6) |  lowbit(dat[i]);
        }
        return N;
    }
    int prev(int val) const {  // (include val)
        if (val < 0) return -1;
        if (val >= N) val = N - 1;
        int k = val >> 6;
        if ((val & 63) < 63) {
            u64 x = dat[k];
            x &= (u64(1) << ((val & 63) + 1)) - 1;
            if (x) return (k << 6) | topbit(x);
            --k;
        }
        for (int i = k; i >= 0; i--) {
            if (dat[i] == 0) continue;
            return (i << 6) | topbit(dat[i]);
        }
        return -1;
    }
    int get_min() const { return next(0); }
    int get_max() const { return prev(N); }

    // input/output operators
    friend istream& operator >> (istream &is, DynamicBitset &db) {
        string S;
        is >> S;
        db.assign(S.size(), 0);
        for (int i = 0; i < S.size(); i++) {
            if (S[i] == '1') db.dat[i >> 6] |= u64(1) << (i & 63);
        }
        return is;
    }
    friend ostream& operator << (ostream &os, const DynamicBitset &db) {
        return os << db.to_string();
    }

    // comparison operators
    bool operator == (const DynamicBitset &r) const {
        assert(size() == r.size());
        for (int i = 0; i < (int)dat.size(); i++) if (dat[i] != r.dat[i]) return false;
        return true;
    }
    bool operator != (const DynamicBitset &r) const {
        assert(size() == r.size());
        return !((*this) == r);
    }

    // set (change to 1)
    void set() {
        fill(dat.begin(), dat.end(), u64(-1));
        resize(N);
    }
    void set(int i) {
        (*this)[i] = 1;
    }
    void set(int L, int R) {
        assert(L <= R);
        while (L < R && (L & 63)) set(L++);
        while (L < R && (R & 63)) set(--R);
        for (int i = (L >> 6); i < (R >> 6); i++) dat[i] = u64(-1);
    }

    // reset (change to 0)
    void reset() {
        fill(dat.begin(), dat.end(), u64(0));
        resize(N);
    }
    void reset(int i) { 
        (*this)[i] = 0;
    }
    void reset(int L, int R) {
        assert(L <= R);
        while (L < R && (L & 63)) reset(L++);
        while (L < R && (R & 63)) reset(--R);
        for (int i = (L >> 6); i < (R >> 6); i++) dat[i] = u64(0);
    }

    // flip (change 0-1)
    void flip() {
        for (int i = 0; i < (int)dat.size() - 1; i++) dat[i] = u64(-1) ^ dat[i];
        int i = (int)dat.size() - 1;
        for (int k = 0; k < 64; k++) {
            if (i * 64 + k >= size()) break;
            flip(i * 64 + k);
        }
    }
    void flip(int i) {
        (*this)[i].flip();
    }
    void flip(int L, int R) {
        assert(L <= R);
        while (L < R && (L & 63)) flip(L++);
        while (L < R && (R & 63)) flip(--R);
        for (int i = (L >> 6); i < (R >> 6); i++) dat[i] ^= u64(-1);
    }

    // logic operators
    DynamicBitset &operator &= (const DynamicBitset &r) {
        assert(size() == r.size());
        for (int i = 0; i < (int)dat.size(); i++) dat[i] &= r.dat[i];
        return *this;
    }
    DynamicBitset &operator |= (const DynamicBitset &r) {
        assert(size() == r.size());
        for (int i = 0; i < (int)dat.size(); i++) dat[i] |= r.dat[i];
        return *this;
    }
    DynamicBitset &operator ^= (const DynamicBitset &r) {
        assert(size() == r.size());
        for (int i = 0; i < (int)dat.size(); i++) dat[i] ^= r.dat[i];
        return *this;
    }
    DynamicBitset operator & (const DynamicBitset &r) const {
        return DynamicBitset(*this) &= r;
    }
    DynamicBitset operator | (const DynamicBitset &r) const {
        return DynamicBitset(*this) |= r;
    }
    DynamicBitset operator ^ (const DynamicBitset &r) const {
        return DynamicBitset(*this) ^= r;
    }
    DynamicBitset operator ~ () const {
        DynamicBitset res(*this);
        res.flip();
        return res;
    }

    // slice [L, R)
    DynamicBitset slice(int L, int R) {
        assert(L <= R);
        DynamicBitset res(R - L);
        int rm = (R - L) & 63;
        for (int i = 0; i < rm; i++) {
            res[R - L - 1] = bool((*this)[R - 1]);
            --R;
        }
        int n = (R - L) >> 6, high = L & 63, s = L >> 6;
        if (!high) {
            for (int i = 0; i < n; i++) {
                res.dat[i] = dat[i + s];
            }
        } else {
            for (int i = 0; i < n; i++) {
                res.dat[i] ^= (dat[i + s] >> high) ^ (dat[i + s + 1] << (64 - high));
            }
        }
        return res;
    }
    
    // apply another DynamicBitset for range [L, R)
    void apply(int L, int R, const DynamicBitset &db) {
        assert(db.size() == R - L);
        int a = 0, b = db.size();
        while (L < R && (L & 63)) (*this)[L++] = bool(db[a++]);
        while (L < R && (R & 63)) (*this)[--R] = bool(db[--b]);
        int l = L >> 6, r = R >> 6, s = a >> 6, high = a & 63;
        if (!high) {
            for (int i = 0; i < r - l; i++) {
                dat[i + l] = db.dat[i + s];
            }
        } else {
            for (int i = 0; i < r - l; i++) {
                dat[i + l] = (db.dat[i + s] >> high) | (db.dat[i + s + 1] << (64 - high));
            }
        }
    }
    void apply_and(int L, int R, const DynamicBitset &db) {
        assert(db.size() == R - L);
        int a = 0, b = db.size();
        while (L < R && (L & 63)) {
            if (!db[a]) (*this)[L] = 0;
            a++, L++;
        }
        while (L < R && (R & 63)) {
            --b, --R;
            if (!db[b]) (*this)[R] = 0;
        }
        int l = L >> 6, r = R >> 6, s = a >> 6, high = a & 63;
        if (!high) {
            for (int i = 0; i < r - l; i++) {
                dat[i + l] &= db.dat[i + s];
            }
        } else {
            for (int i = 0; i < r - l; i++) {
                dat[i + l] &= (db.dat[i + s] >> high) | (db.dat[i + s + 1] << (64 - high));
            }
        }
    }
    void apply_or(int L, int R, const DynamicBitset &db) {
        assert(db.size() == R - L);
        int a = 0, b = db.size();
        while (L < R && (L & 63)) {
            dat[L >> 6] |= u64(db[a]) << (L & 63);
            ++a, ++L;
        }
        while (L < R && (R & 63)) {
            --b, --R;
            dat[R >> 6] |= u64(db[b]) << (R & 63);
        }
        int l = L >> 6, r = R >> 6, s = a >> 6, high = a & 63;
        if (!high) {
            for (int i = 0; i < r - l; i++) {
                dat[i + l] |= db.dat[i + s];
            }
        } else {
            for (int i = 0; i < r - l; i++) {
                dat[i + l] |= (db.dat[i + s] >> high) | (db.dat[i + s + 1] << (64 - high));
            }
        }
    }
    void apply_xor(int L, int R, const DynamicBitset &db) {
        assert(db.size() == R - L);
        int a = 0, b = db.size();
        while (L < R && (L & 63)) {
            dat[L >> 6] ^= u64(db[a]) << (L & 63);
            ++a, ++L;
        }
        while (L < R && (R & 63)) {
            --b, --R;
            dat[R >> 6] ^= u64(db[b]) << (R & 63);
        }
        int l = L >> 6, r = R >> 6, s = a >> 6, high = a & 63;
        if (!high) {
            for (int i = 0; i < r - l; i++) {
                dat[i + l] ^= db.dat[i + s];
            }
        } else {
            for (int i = 0; i < r - l; i++) {
                dat[i + l] ^= (db.dat[i + s] >> high) | (db.dat[i + s + 1] << (64 - high));
            }
        }
    }

    // apply xor for [L, N) (use in matrix sweeping)
    void apply_xor(int L, const DynamicBitset &db) {
        assert(size() == db.size());
        assert(0 <= L && L < size());
        for (int i = (L >> 6); i < (int)dat.size(); i++) dat[i] ^= db.dat[i];
    }
};
string DynamicBitset::CONV[256];


//------------------------------//
// Mod 2 Matrix
//------------------------------//

// binary matrix
struct BinaryMatrix {
    // inner value
    int H, W;
    vector<DynamicBitset> val;

    // constructors
    BinaryMatrix() {}
    BinaryMatrix(const BinaryMatrix&) = default;
    BinaryMatrix& operator = (const BinaryMatrix&) = default;
    BinaryMatrix(int h, int w, int x = 0) : H(h), W(w), val(h, DynamicBitset(w, x)) {}
    void init(int h, int w, int x = 0) {
        H = h, W = w;
        val.assign(h, DynamicBitset(w, x));
    }
    void resize(int h, int w) {
        H = h, W = w;
        val.resize(h);
        for (int i = 0; i < h; ++i) val[i].resize(w);
    }

    // getter and debugger
    constexpr int height() const { return H; }
    constexpr int width() const { return W; }
    constexpr bool empty() const { return height() == 0; }
    DynamicBitset& operator [] (int i) { return val[i]; }
    const DynamicBitset& operator [] (int i) const { return val[i]; }
    friend ostream& operator << (ostream &os, const BinaryMatrix &mat) {
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) {
                if (j) os << ' ';
                os << mat.val[i][j];
            }
            os << '\n';
        }
        return os;
    }
    
    // comparison operators
    bool operator == (const BinaryMatrix &r) const {
        return this->val == r.val;
    }
    bool operator != (const BinaryMatrix &r) const {
        return this->val != r.val;
    }

    // arithmetic operators
    const BinaryMatrix& operator += (const BinaryMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i) val[i] ^= r.val[i];
        return *this;
    }
    const BinaryMatrix& operator -= (const BinaryMatrix &r) {
        assert(height() == r.height());
        assert(width() == r.width());
        for (int i = 0; i < height(); ++i) val[i] ^= r.val[i];
        return *this;
    }
    const BinaryMatrix& operator *= (int v) {
        assert(v == 0 || v == 1);
        if (v == 1) return *this;
        else for (int i = 0; i < height(); ++i) val[i].reset();
        return *this;
    }
    const BinaryMatrix& operator *= (const BinaryMatrix &r) {
        assert(width() == r.height());
        BinaryMatrix res(height(), r.width());
        BinaryMatrix tr = r.trans();
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < r.width(); ++j)
                res[i][j] = val[i].dot_mod2(tr[j]);
        return (*this) = res;
    }
    BinaryMatrix operator + () const { 
        return BinaryMatrix(*this);
    }
    BinaryMatrix operator - () const {
        return BinaryMatrix(*this);
    }
    BinaryMatrix operator + (const BinaryMatrix &r) const { 
        return BinaryMatrix(*this) += r;
    }
    BinaryMatrix operator - (const BinaryMatrix &r) const {
        return BinaryMatrix(*this) -= r;
    }
    BinaryMatrix operator * (int v) const {
        return BinaryMatrix(*this) *= v;
    }
    BinaryMatrix operator * (const BinaryMatrix &r) const {
        return BinaryMatrix(*this) *= r;
    }
    DynamicBitset operator * (const DynamicBitset &v) const {
        assert(width() == v.size());
        DynamicBitset res(height(), 0);
        for (int i = 0; i < height(); i++) res[i] = val[i].dot_mod2(v);
        return res;
    }

    // transpose
    BinaryMatrix trans() const {
        BinaryMatrix res(width(), height());
        for (int row = 0; row < width(); row++)
            for (int col = 0; col < height(); col++)
                res[row][col] = val[col][row];
        return res;
    }
    friend BinaryMatrix trans(const BinaryMatrix &mat) {
        return mat.trans();
    }
    
    // pow
    BinaryMatrix pow(long long n) const {
        assert(height() == width());
        BinaryMatrix res(height(), width());
        BinaryMatrix mul(*this);
        for (int row = 0; row < height(); ++row) res[row][row] = 1;
        while (n > 0) {
            if (n & 1) res = res * mul;
            mul = mul * mul;
            n >>= 1;
        }
        return res;
    }
    friend BinaryMatrix pow(const BinaryMatrix &mat, long long n) {
        return mat.pow(n);
    }

    // gauss-jordan
    int find_pivot(int cur_rank, int col) const {
        int pivot = -1;
        for (int row = cur_rank; row < height(); ++row) {
            if (val[row][col] != 0) {
                pivot = row;
                break;
            }
        }
        return pivot;
    }
    void sweep(int cur_rank, int col, int pivot, bool sweep_upper = true) {
        swap(val[pivot], val[cur_rank]);
        assert(val[cur_rank][col] != 0);
        int row_start = (sweep_upper ? 0 : cur_rank + 1);
        for (int row = row_start; row < height(); ++row) {
            if (row != cur_rank && val[row][col] != 0) {
                auto fac = val[row][col];
                val[row].apply_xor(cur_rank, val[cur_rank]);
            }
        }
    }
    int gauss_jordan(int not_sweep_width = 0, bool sweep_upper = true) {
        int rank = 0;
        for (int col = 0; col < width(); ++col) {
            if (col == width() - not_sweep_width) break;
            int pivot = find_pivot(rank, col);
            if (pivot == -1) continue;
            sweep(rank++, col, pivot, sweep_upper);
        }
        return rank;
    }
    int gauss_jordan(vector<int> &core, int not_sweep_width, bool sweep_upper = true) {
        core.clear();
        int rank = 0;
        for (int col = 0; col < width(); ++col) {
            if (col == width() - not_sweep_width) break;
            int pivot = find_pivot(rank, col);
            if (pivot == -1) continue;
            core.push_back(col);
            sweep(rank++, col, pivot, sweep_upper);
        }
        return rank;
    }
    friend int gauss_jordan(BinaryMatrix &mat, int not_sweep_width = 0, bool sweep_upper = true) {
        return mat.gauss_jordan(not_sweep_width, sweep_upper);
    }

    // rank
    int get_rank() const {
        if (height() == 0 || width() == 0) return 0;
        BinaryMatrix A(*this);
        return A.gauss_jordan(0, false);
    }
    friend int get_rank(const BinaryMatrix &mat) {
        return mat.get_rank();
    }

    // find one solution
    friend int linear_equation
    (const BinaryMatrix &mat, const DynamicBitset &b, DynamicBitset &res) {
        assert(mat.height() == b.size());

        // extend
        BinaryMatrix A(mat.height(), mat.width() + 1);
        for (int i = 0; i < mat.height(); ++i) {
            A[i].apply(0, mat.width(), mat[i]);
            A[i][A[i].size()-1] = b[i];
        }
        int rank = A.gauss_jordan(1);
        
        // check if it has no solution
        for (int row = rank; row < mat.height(); ++row) {
            if (A[row].back() != 0) return -1;
        }

        // answer
        res.assign(mat.width(), 0);
        for (int i = 0; i < rank; ++i) res[i] = A[i].back();
        return rank;
    }
    friend int linear_equation(const BinaryMatrix &mat, const DynamicBitset &b) {
        DynamicBitset res;
        return linear_equation(mat, b, res);
    }

    // find all solutions
    friend int linear_equation
    (const BinaryMatrix &mat, const DynamicBitset &b, DynamicBitset &res, vector<DynamicBitset> &zeros) {
        // extend
        BinaryMatrix A(mat.height(), mat.width() + 1);
        for (int i = 0; i < mat.height(); ++i) {
            A[i].apply(0, mat.width(), mat[i]);
            A[i][A[i].size()-1] = b[i];
        }
        vector<int> core;
        int rank = A.gauss_jordan(core, 1);
        
        // check if it has no solution
        for (int row = rank; row < mat.height(); ++row) {
            if (A[row].back() != 0) return -1;
        }

        // construct the core solution
        res.assign(mat.width(), 0);
        for (int i = 0; i < (int)core.size(); i++) res[core[i]] = A[i].back();
    
        // construct the all solutions
        zeros.clear();
        DynamicBitset used(mat.width(), 0);
        for (auto c : core) used[c] = 1;
        for (int j = 0; j < mat.width(); j++) {
            if (used[j]) continue;
            DynamicBitset zero(mat.width(), 0);
            zero[j] = 1;
            for (int i = 0; i < (int)core.size(); i++) zero[core[i]] = -A[i][j];
            zeros.push_back(zero);
        }
        return rank;
    }
    
    // determinant
    int det() const {
        assert(height() == width());
        if (height() == 0) return 1;
        BinaryMatrix A(*this);
        int rank = 0, res = 1;
        for (int col = 0; col < width(); ++col) {
            int pivot = A.find_pivot(rank, col);
            if (pivot == -1) return 0;
            res *= A[pivot][rank];
            A.sweep(rank++, col, pivot, false);
        }
        return res;
    }
    friend int det(const BinaryMatrix &mat) {
        return mat.det();
    }

    // inv
    BinaryMatrix inv() const {
        assert(height() == width());

        // extend
        BinaryMatrix A(height(), width() + height());
        for (int i = 0; i < height(); ++i) {
            A[i].apply(0, width(), val[i]);
            A[i][i + width()] = 1;
        }
        vector<int> core;
        int rank = A.gauss_jordan(height(), true);

        // gauss jordan
        if (rank < height()) return BinaryMatrix(0, 0);
        BinaryMatrix res(height(), width());
        for (int i = 0; i < height(); ++i) {
            res[i].apply(0, width(), A[i].slice(width(), width() + height()));
        }
        return res;
    }
    friend BinaryMatrix inv(const BinaryMatrix &mat) {
        return mat.inv();
    }
};


//------------------------------//
// Exapmles
//------------------------------//

// Yosupo Library Checker - Matrix Product (Mod 2)
void Yosupo_Matrix_Product_Mod2() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N, M, K;
    cin >> N >> M >> K;
    BinaryMatrix A(N, M), B(M, K);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < M; i++) cin >> B[i];
    auto C = A * B;
    for (int i = 0; i < C.height(); i++) cout << C[i] << '\n';
}

// Yosupo Library Checker - Determinant of Matrix (Mod 2)
void Yosupo_Determinant_of_Matrix_Mod2() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N;
    cin >> N;
    BinaryMatrix A(N, N);
    for (int i = 0; i < N; i++) cin >> A[i];
    auto res = det(A);
    cout << res << '\n';
}

// Yosupo Library Checker - Rank of Matrix (Mod 2)
void Yosupo_Rank_of_Matrix_Mod2() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N, M;
    cin >> N >> M;
    if (N == 0 || M == 0) {
        cout << 0 << endl;
        return;
    }
    bool flip = false;
    if (N > M) swap(N, M), flip = true;
    BinaryMatrix A(N, M);
    if (!flip) {
        for (int i = 0; i < N; i++) cin >> A[i];
    } else {
        for (int i = 0; i < M; i++) {
            string S;
            cin >> S;
            for (int j = 0; j < N; j++) {
                if (S[j] == '1') A[j].set(i);
            }
        }
    }
    auto res = get_rank(A);
    cout << res << '\n';
}

// Yosupo Library Checker - System of Linear Equations (Mod 2)
void Yosupo_System_of_Linear_Equations_Mod2() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N, M;
    cin >> N >> M;
    BinaryMatrix A(N, M);
    for (int i = 0; i < N; i++) cin >> A[i];
    DynamicBitset b(N), res;
    cin >> b;
    vector<DynamicBitset> zeros;
    int rank = linear_equation(A, b, res, zeros);
    if (rank == -1) {
        cout << -1 << endl;
        return;
    }
    cout << zeros.size() << endl;
    cout << res << endl;
    for (int i = 0; i < (int)zeros.size(); i++) cout << zeros[i] << endl;
}

// Yosupo Library Checker - Inverse Matrix (Mod 2)
void Yosupo_Inverse_Matrix_Mod2() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N;
    cin >> N;
    BinaryMatrix A(N, N);
    for (int i = 0; i < N; i++) cin >> A[i];
    auto res = inv(A);
    if (res.empty()) {
        cout << -1 << endl;
        return;
    }
    for (int i = 0; i < (int)res.height(); i++) cout << res[i] << endl;
}


int main() {
    //Yosupo_Matrix_Product_Mod2();
    //Yosupo_Determinant_of_Matrix_Mod2();
    Yosupo_Rank_of_Matrix_Mod2();
    //Yosupo_System_of_Linear_Equations_Mod2();
    //Yosupo_Inverse_Matrix_Mod2();
}