//
// 動的 bitset
//   ・可変長にしたいとき
//   ・特定の範囲を取り出して操作したいとき
//
// verified:
//   Yosupo Library Checker - Matrix Product (Mod 2)
//     https://judge.yosupo.jp/problem/matrix_product_mod_2
//
//   Codeforces Round 1079 (Div. 1) E2. Fuzzy Concatenation (Hard version)
//     https://codeforces.com/contest/2196/problem/E2
//


#include <bits/stdc++.h>
using namespace std;


//------------------------------//
// Utility
//------------------------------//

template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl

// num of i such that (x & (1 << i)) != 0
int popcnt(int x) { return __builtin_popcount(x); }
int popcnt(unsigned int x) { return __builtin_popcount(x); }
int popcnt(long long x) { return __builtin_popcountll(x); }
int popcnt(unsigned long long x) { return __builtin_popcountll(x); }
int popcnt_mod2(int x) { return __builtin_parity(x); }
int popcnt_mod2(unsigned int x) { return __builtin_parity(x); }
int popcnt_mod2(long long x) { return __builtin_parityll(x); }
int popcnt_mod2(unsigned long long x) { return __builtin_parityll(x); }


//------------------------------//
// 動的 bitset
//------------------------------//

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
// Examples
//------------------------------//

void small_test() {
    DynamicBitset a(20), b("00000000001111000010");
    a.set(3), a.set(8, 11);
    assert(a == DynamicBitset("00010000111000000000"));
    assert(b == DynamicBitset("00000000001111000010"));
    assert((a & b) == DynamicBitset("00000000001000000000"));
    assert((a | b) == DynamicBitset("00010000111111000010"));
    assert((a ^ b) == DynamicBitset("00010000110111000010"));
    a[1] = 1;
    assert(a == DynamicBitset("01010000111000000000"));
    vector<int> prev{-1, 1, 1, 3, 3, 3, 3, 3, 8, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
    vector<int> next{1, 1, 3, 3, 8, 8, 8, 8, 8, 9, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20};
    for (int i = 0; i < 20; i++) {
        assert(a.prev(i) == prev[i]);
        assert(a.next(i) == next[i]);
    }
    a.flip(10, 17);
    assert(a == DynamicBitset("01010000110111111000"));
    assert(a.count() == 10);
    assert(a.count(3, 16) == 8);
    auto c = a.slice(6, 19);
    c[7] = 0;
    assert(c == DynamicBitset("0011011011100"));
    a.apply(1, 1+c.size(), c);
    assert(a == DynamicBitset("00011011011100111000"));
    a.apply_and(2, 2+c.size(), c);
    assert(a == DynamicBitset("00001001001100011000"));
    a.apply_or(0, c.size(), c);
    assert(a == DynamicBitset("00111111111100011000"));
    a.apply_xor(20-c.size(), 20, c);
    assert(a == DynamicBitset("00111111100111000100"));
    a[5] = b[3];
    assert(a == DynamicBitset("00111011100111000100"));
}


int main() {
    small_test();
}